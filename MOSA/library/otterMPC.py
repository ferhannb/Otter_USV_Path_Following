#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
otter.py: 
    Class for the Maritime Robotics Otter USV, www.maritimerobotics.com. 
    The length of the USV is L = 2.0 m. The constructors are:

    otter()                                          
        Step inputs for propeller revolutions n1 and n2
    otter('headingAutopilot',psi_d,V_current,beta_current,tau_X)  
       Heading autopilot with options:
          psi_d: desired yaw angle (deg)
          V_current: current speed (m/s)
          beta_c: current direction (deg)
          tau_X: surge force, pilot input (N)
        
Methods:
    
[nu,u_actual] = dynamics(eta,nu,u_actual,u_control,sampleTime) returns 
    nu[k+1] and u_actual[k+1] using Euler's method. The control inputs are:

    u_control = [ n1 n2 ]' where 
      n1: propeller shaft speed, left (rad/s)
      n2: propeller shaft speed, right (rad/s)

u = headingAutopilot(eta,nu,sampleTime) 
    PID controller for automatic heading control based on pole placement.

u = stepInput(t) generates propeller step inputs.

[n1, n2] = controlAllocation(tau_X, tau_N)     
    Control allocation algorithm.
    
References: 
  T. I. Fossen (2021). Handbook of Marine Craft Hydrodynamics and Motion 
     Control. 2nd. Edition, Wiley. 
     URL: www.fossen.biz/wiley            

Author:     Thor I. Fossen
"""


import numpy as np
import math
import time


#from OtterSimulator.control import PIDpolePlacement
#from OtterSimulator.gnc import Smtrx, Hmtrx, m2c, crossFlowDrag, sat, attitudeEuler
from gnc import Smtrx, Hmtrx, m2c, crossFlowDrag, sat, attitudeEuler
from path_track_2 import LineofSight
from R_calculation import R_Calculator
from otterradius import speed_generate
import time
# from mpccontroller import MPCUSV
import casadi as ca
# Class Vehicle
class Otter:
    """
    otter()                          Step inputs for n1 and n2
    otter('headingAutopilot',psi_d)  Heading autopilot, desired yaw angle (deg)
    """

    def __init__(self, 
                controlSystem="stepInput", 
                r = 0, 
                V_current = 0, 
                beta_current = 0, 
                tau_X = 120):

        if controlSystem == "headingAutopilot":
            self.controlDescription = (
                "Heading autopilot, psi_d = "
                + str(r)
                + " deg"
                )
        else:
            self.controlDescription = "Step inputs for n1 and n2"
            controlSystem = "stepInput"
        self.C = np.zeros((6,6))
        self.ref = r
        self.V_c = V_current
        self.beta_c = beta_current
        self.controlMode = controlSystem
        self.tauX = tau_X  # surge force (N)

        # Initialize the Otter USV model
        self.T_n = 0.32  # propeller time constants (s)
        self.L = 2.0    # Length (m)
        self.B = 1.08   # beam (m)
        self.nu = np.array([0, 0, 0, 0, 0, 0], float)  # velocity vector
        self.u_actual = np.array([0, 0], float)  # propeller revolution states
        self.name = "Otter USV (see 'otter.py' for more details)"

        self.controls = [
            "Left propeller shaft speed (rad/s)",
            "Right propeller shaft speed (rad/s)",
        ]
        self.dimU = len(self.controls)

        # Constants
        g = 9.81    # acceleration of gravity (m/s^2)
        rho = 1025  # density of water

        m = 55.0    # mass (kg)
        mp = 25.0   # Payload (kg)
        self.m_total = m + mp
        rp = np.array([0, 0, 0], float)     # location of payload (m)
        rg = np.array([0.2, 0, -0.2], float)    # CG for hull only (m)
        rg = (m * rg + mp * rp) / (m + mp)      # CG corrected for payload
        self.S_rg = Smtrx(rg)
        self.H_rg = Hmtrx(rg)
        self.S_rp = Smtrx(rp)

        R44 = 0.4 * self.B  # radii of gyration (m)
        R55 = 0.25 * self.L
        R66 = 0.25 * self.L
        T_yaw = 1.0         # time constant in yaw (s)
        Umax = 6 * 0.5144   # max forward speed (m/s)

        # Data for one pontoon
        self.B_pont = 0.25  # beam of one pontoon (m)
        y_pont = 0.395      # distance from centerline to waterline centroid (m)
        Cw_pont = 0.75      # waterline area coefficient (-)
        Cb_pont = 0.4       # block coefficient, computed from m = 55 kg

        # Inertia dyadic, volume displacement and draft
        nabla = (m + mp) / rho  # volume
        self.T = nabla / (2 * Cb_pont * self.B_pont * self.L)  # draft
        Ig_CG = m * np.diag(np.array([R44 ** 2, R55 ** 2, R66 ** 2]))
        self.Ig = Ig_CG - m * self.S_rg @ self.S_rg - mp * self.S_rp @ self.S_rp

        # Experimental propeller data including lever arms
        self.l1 = -y_pont  # lever arm, left propeller (m)
        self.l2 = y_pont  # lever arm, right propeller (m)
        self.k_pos = 0.02216 / 2  # Positive Bollard, one propeller
        self.k_neg = 0.01289 / 2  # Negative Bollard, one propeller
        self.n_max = math.sqrt((0.5 * 24.4 * g) / self.k_pos)  # max. prop. rev.
        self.n_min = -math.sqrt((0.5 * 13.6 * g) / self.k_neg) # min. prop. rev.

        # MRB_CG = [ (m+mp) * I3  O3      (Fossen 2021, Chapter 3)
        #               O3       Ig ]
        MRB_CG = np.zeros((6, 6))
        MRB_CG[0:3, 0:3] = (m + mp) * np.identity(3)
        MRB_CG[3:6, 3:6] = self.Ig
        MRB = self.H_rg.T @ MRB_CG @ self.H_rg

        # Hydrodynamic added mass (best practice)
        Xudot = -0.1 * m
        Yvdot = -1.5 * m
        Zwdot = -1.0 * m
        Kpdot = -0.2 * self.Ig[0, 0]
        Mqdot = -0.8 * self.Ig[1, 1]
        Nrdot = -1.7 * self.Ig[2, 2]

        self.MA = -np.diag([Xudot, Yvdot, Zwdot, Kpdot, Mqdot, Nrdot])

        # System mass matrix
        self.M = MRB + self.MA
        
        self.Minv = np.linalg.inv(self.M)

        # Hydrostatic quantities (Fossen 2021, Chapter 4)
        Aw_pont = Cw_pont * self.L * self.B_pont  # waterline area, one pontoon
        I_T = (
            2
            * (1 / 12)
            * self.L
            * self.B_pont ** 3
            * (6 * Cw_pont ** 3 / ((1 + Cw_pont) * (1 + 2 * Cw_pont)))
            + 2 * Aw_pont * y_pont ** 2
        )
        I_L = 0.8 * 2 * (1 / 12) * self.B_pont * self.L ** 3
        KB = (1 / 3) * (5 * self.T / 2 - 0.5 * nabla / (self.L * self.B_pont))
        BM_T = I_T / nabla  # BM values
        BM_L = I_L / nabla
        KM_T = KB + BM_T    # KM values
        KM_L = KB + BM_L
        KG = self.T - rg[2]
        GM_T = KM_T - KG    # GM values
        GM_L = KM_L - KG

        G33 = rho * g * (2 * Aw_pont)  # spring stiffness
        G44 = rho * g * nabla * GM_T
        G55 = rho * g * nabla * GM_L
        G_CF = np.diag([0, 0, G33, G44, G55, 0])  # spring stiff. matrix in CF
        LCF = -0.2
        H = Hmtrx(np.array([LCF, 0.0, 0.0]))  # transform G_CF from CF to CO
        self.G = H.T @ G_CF @ H

        # Natural frequencies
        w3 = math.sqrt(G33 / self.M[2, 2])
        w4 = math.sqrt(G44 / self.M[3, 3])
        w5 = math.sqrt(G55 / self.M[4, 4])

        # Linear damping terms (hydrodynamic derivatives)
        Xu = -24.4 * g / Umax  # specified using the maximum speed
        Yv = 0
        Zw = -2 * 0.3 * w3 * self.M[2, 2]  # specified using relative damping
        Kp = -2 * 0.2 * w4 * self.M[3, 3]
        Mq = -2 * 0.4 * w5 * self.M[4, 4]
        Nr = -self.M[5, 5] / T_yaw  # specified by the time constant T_yaw

        self.D = -np.diag([Xu, Yv, Zw, Kp, Mq, Nr])

        # Trim: theta = -7.5 deg corresponds to 13.5 cm less height aft
        self.trim_moment = 0
        self.trim_setpoint = 280

        # Propeller configuration/input matrix
        B = self.k_pos * np.array([[1, 1], [-self.l1, -self.l2]])
        self.Binv = np.linalg.inv(B)

        # Heading autopilot
        self.e_int = 0  # integral state
        self.wn = 1.2  # PID pole placement
        self.zeta = 0.8

        # Reference model
        self.r_max = 10 * math.pi / 180  # maximum yaw rate
        self.psi_d = 0  # angle, angular rate and angular acc. states
        self.r_d = 0
        self.a_d = 0
        self.wn_d = self.wn / 5  # desired natural frequency in yaw
        self.zeta_d = 1  # desired relative damping ratio


    def dynamics(self, eta, nu, u_actual, u_control, sampleTime):
        """
        [nu,u_actual] = dynamics(eta,nu,u_actual,u_control,sampleTime) integrates
        the Otter USV equations of motion using Euler's method.
        """

        # Input vector
        n = np.array([u_actual[0], u_actual[1]])
        #print('n:',n)
        # Current velocities
        u_c = self.V_c * math.cos(self.beta_c - eta[5])  # current surge vel.
        v_c = self.V_c * math.sin(self.beta_c - eta[5])  # current sway vel.

        nu_c = np.array([u_c, v_c, 0, 0, 0, 0], float)  # current velocity vector
        nu_r = nu - nu_c  # relative velocity vector

        # Rigid body and added mass Coriolis and centripetal matrices
        # CRB_CG = [ (m+mp) * Smtrx(nu2)          O3   (Fossen 2021, Chapter 6)
        #              O3                   -Smtrx(Ig*nu2)  ]
        CRB_CG = np.zeros((6, 6))
        CRB_CG[0:3, 0:3] = self.m_total * Smtrx(nu[3:6])
        CRB_CG[3:6, 3:6] = -Smtrx(np.matmul(self.Ig, nu[3:6]))
        CRB = self.H_rg.T @ CRB_CG @ self.H_rg  # transform CRB from CG to CO
        # print('Crb',CRB_CG)
        CA = m2c(self.MA, nu_r)
        CA[5, 0] = 0  # assume that the Munk moment in yaw can be neglected
        CA[5, 1] = 0  # if nonzero, must be balanced by adding nonlinear damping

        self.C = CRB + CA

        # Ballast
        g_0 = np.array([0.0, 0.0, 0.0, 0.0, self.trim_moment, 0.0])

        # Control forces and moments - with propeller revolution saturation
        thrust = np.zeros(2)
        for i in range(0, 2):

            n[i] = sat(n[i], self.n_min, self.n_max)  # saturation, physical limits

            if n[i] > 0:  # positive thrust
                thrust[i] = self.k_pos * n[i] * abs(n[i])
            else:  # negative thrust
                thrust[i] = self.k_neg * n[i] * abs(n[i])

        # Control forces and moments
        tau = np.array(
            [
                thrust[0] + thrust[1],
                0,
                0,
                0,
                0,
                -self.l1 * thrust[0] - self.l2 * thrust[1],
            ]
        )
        

        # Hydrodynamic linear damping + nonlinear yaw damping
        tau_damp = -np.matmul(self.D, nu_r)
        tau_damp[5] = tau_damp[5] - 10 * self.D[5, 5] * abs(nu_r[5]) * nu_r[5]
        # print('tau_damp',tau_damp)
        # State derivatives (with dimension)
        tau_crossflow = crossFlowDrag(self.L, self.B_pont, self.T, nu_r)
        sum_tau = (
            tau
            + tau_damp
            + tau_crossflow
            - np.matmul(self.C, nu_r)
            - np.matmul(self.G, eta)
            - g_0
        )


        nu_dot = np.matmul(self.Minv, sum_tau)  # USV dynamics

        n_dot = (u_control - n) / self.T_n  # propeller dynamics
        trim_dot = (self.trim_setpoint - self.trim_moment) / 5  # trim dynamics

        # Forward Euler integration [k+1]
        nu = nu + sampleTime * nu_dot
        n = n + sampleTime * n_dot
        self.trim_moment = self.trim_moment + sampleTime * trim_dot

        u_actual = np.array(n, float)
       

        

        ##### MPC 3 DOF MODEL #########

        # self.M_inv_3 = self.Minv[[0,1,5],:][:,[0,1,5]]
        # tau_3 = np.transpose(tau[[0,1,5]])
        # tau_damp_3 = np.transpose(tau_damp[[0,1,5]])
        # tau_crossflow_3 = np.transpose(tau_crossflow[[0,1,5]])
        # self.D_3 = self.D[[0,1,5],:][:,[0,1,5]]
        # self.C_3 = self.C[[0,1,5],:][:,[0,1,5]] 
        # G_3 = self.G[[0,1,5],:][:,[0,1,5]]
        # g0_3 = np.transpose([g_0[[0,1,5]]])
        # nu_r_3 = np.transpose(nu_r[[0,1,5]])
        # eta_3 =np.transpose(eta[[0,1,5]])
        
        
        

        return nu, u_actual

    def dynamic_eq(self):
        self.M_inv_3 = self.Minv#[[0,1,5],:][:,[0,1,5]]
        self.D_3 = self.D#[[0,1,5],:][:,[0,1,5]]
        self.C_3 = self.C#[[0,1,5],:][:,[0,1,5]]

        return self.D_3,self.C_3,self.M_inv_3

       
    
    def initialize_otter(self, sample_time):
        self.sample_time=sample_time
        self.nu = np.array([0, 0, 0, 0, 0, 0], float)
        self.u_actual = np.array([0, 0], float)
        self.u_control = np.array([0, 0], float)
        self.current_eta = np.array([0, 0, 0, 0, 0, 0], float)

        return
  

    def stepInput(self, iteration, length):   #vehicle inputs 
        """
        u = stepInput(t) generates propeller step inputs.
        """
        n1 = 20.94# rad/s
        n2 = -20.94# rad/s

        # if iteration > int(length/2):
        #     n1 = 50
        #     n2 = 50

        u_control = np.array([n1, n2], float)

        return u_control

##################################################################################################
"""sampling time 0.1'den büyük olmamalı (default=0.02)
    input değerleri
    -------------------------
    current_eta=6 DOF konumlar
    nu=6 DOF hizlar
    u_actual=pervane hizi rad/s pervane (modelinin çalışması için)
    u_control=pervane hiz referans değeri rad/s
    ****************
    output değerleri
    -------------------------------------------
    current_eta: X-Y de hız grafilerini çizdirebilmek için
    heading : anlık dünya koordinatına göre heading 
    velocity: aracın gittiği u-v bileşke hızı
    u-actual: pervane modelinin çalışması için sisteme tekrar girecek pervane rpm değeri
    
    """

# def start(steptime):
#     vehicle=Otter()
#     global current_eta
#     global nu
#     global u_actual
#     global u_control
#     global sampleTime
#     global speed
#     global heading

    
#     nu = np.array([0, 0, 0, 0, 0, 0], float)
#     u_actual = np.array([0, 0], float)
#     u_control = np.array([0, 0], float)           
#     current_eta = np.array([0, 10, 0, 0, 0, math.radians(0)], float)
#     heading=current_eta[5]*180/math.pi %360
#     speed = math.sqrt(nu[0]**2+nu[1]**2)
#     sampleTime=steptime

#     return

# def function(u_control):
#     vehicle=Otter()
#     global current_eta
#     global nu
#     global u_actual
#     global sampleTime
#     global heading
#     global speed

#     sampleTime=0.02
#     [nu, u_actual]=vehicle.dynamics(current_eta, nu, u_actual, u_control, sampleTime)
#     current_eta = attitudeEuler(current_eta, nu, sampleTime)
#     heading=current_eta[5]*180/math.pi %360
#     speed=math.sqrt(nu[0]**2+nu[1]**2)
#     # print('--------------------')

#     # output=[current_eta[0],current_eta[1],u_actual[0],u_actual[1],heading,nu]
#     output = [current_eta,nu,u_actual]
#     time.sleep(sampleTime) # simulasyon süresi ayarlama
#     return output 

# def control_allocation(u_avr,u_diff):
#     max_frw_rpm = 100
#     max_rvs_rpm = -100
#     global u_control

#     if u_avr>max_frw_rpm:
#         u_avr=max_frw_rpm
#     elif u_avr<max_rvs_rpm:
#         u_avr=max_rvs_rpm
#     # print('u_avr---',u_avr)
#     n1=u_avr-u_diff
#     n2=u_avr+u_diff

#     if n1>max_frw_rpm:
#         n1=max_frw_rpm
#     if n1<max_rvs_rpm:
#         n1=max_rvs_rpm
#     if n2>max_frw_rpm:
#         n2=max_frw_rpm   
#     if n2<max_rvs_rpm:
#         n2=max_rvs_rpm

#     u_control=[n1,n2]


# def speed_control(set_point):

#     global nu
#     vehicle_velocity=math.sqrt(nu[0]**2+nu[1]**2)
#     error=set_point-vehicle_velocity
#     if set_point==0:
#             set_point=0.0001
#     else:
#             pass
#     U_f = (3.684*set_point**3-23.44*set_point**2+67.35*set_point+12.3)
#     U_p = 10*error
#     u_avg = U_f+U_p
#     # print('u_avg',u_avg)
#     return u_avg
    



# if __name__ == "__main__": 

#     start(0.02)
#     Wpx = [20,60,100,120]
#     Wpy = [20,80,30,180]
#     los2=LineofSight()
#     R_cal=R_Calculator()

#     mpc = MPCUSV()

#     x_list=[]
#     y_list=[]
#     d_list=[]
#     init_heading=0.0
#     r_list=[]

#     prev_x=0
#     prev_y=0
    
#     prev_heading=0
#     ds = 0.01 #[m]
#     coeff = los2.path_generate(Wpx,Wpy)
#     v_list=[]
#     head_list=[]
#     init_states = [0,0,0,0,10,0]
#     for i in range(4000):
        
#         los_output=los2.execute(current_eta,Wpx,Wpy)
#         mpc_output= mpc.execute_MPC(init_states,[0,0,0,los_output['x_los'],los_output['y_los'],los_output['chi_d']])
  
#         u_d = np.array(mpc_output['control_signal'].full())
#         u_control = [float(u_d[0][0]),float(u_d[1][0])]


#         v_list.append(np.sqrt(nu[0]**2+nu[1]**2))
#         head_list.append(current_eta[5]*180/math.pi %360)
#         x_list.append(current_eta[0])
#         y_list.append(current_eta[1])

#         output=function(u_control)
#         x_init=output[0][0] 
#         y_init=output[0][1]
#         yaw_init=output[0][5]

#         init_states = ca.vertcat(0,0,0,x_init,y_init,yaw_init)
#         radius=R_cal.R_cal(current_eta)
#         los2.los_simulation(current_eta,Wpx,Wpy,u_control)
   
