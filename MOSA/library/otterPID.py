#!/usr/bin/env python
# -*- coding: utf-8 -*-


## Otter PID kontrol tasarımı
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


from calendar import SATURDAY
from stat import UF_IMMUTABLE
from timeit import repeat
import numpy as np
import math
import time
from path_track_3 import LineofSightfirst


#from OtterSimulator.control import PIDpolePlacement
#from OtterSimulator.gnc import Smtrx, Hmtrx, m2c, crossFlowDrag, sat, attitudeEuler
from gnc import Smtrx, Hmtrx, m2c, crossFlowDrag, sat, attitudeEuler
from path_track_2 import LineofSight
from R_calculation import R_Calculator
from otterradius import speed_generate
from path_generation.CubicCurve import generate_curve
from stanley import Stanley





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

        return nu, u_actual


    
    def initialize_otter(self, sample_time):
        self.sample_time=sample_time
        self.nu = np.array([0, 0, 0, 0, 0, 0], float)
        self.u_actual = np.array([0, 0], float)
        self.u_control = np.array([0, 0], float)
        self.current_eta = np.array([10, 20, 0, 0, 0, 0], float)

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
def start(steptime):
    vehicle=Otter()
    global current_eta
    global nu
    global u_actual
    global u_control
    global sampleTime
    global speed
    global heading
    global pre_error_heading
    global pre_error_speed 
    global U_i
    global saturation_limit_speed_max
    global saturation_limit_speed_min
    global U_sat_old

  

    pre_error_heading =0
    pre_error_speed = 0
    U_i = 0
    saturation_limit_speed_max =104
    saturation_limit_speed_min =-104
    U_sat_old = 0

    nu = np.array([0, 0, 0, 0, 0, 0], float)
    u_actual = np.array([0, 0], float)
    u_control = np.array([0, 0], float)           
    current_eta = np.array([10, 10, 0, 0, 0, math.radians(75)], float)
    heading=current_eta[5]*180/math.pi %360
    speed = math.sqrt(nu[0]**2+nu[1]**2)
    sampleTime=steptime

    return

def function(u_control):
    vehicle=Otter()
    global current_eta
    global nu
    global u_actual
    global sampleTime
    global heading
    global speed
    
    sampleTime=0.02
    ttime1=time.time()

    [nu, u_actual]=vehicle.dynamics(current_eta, nu, u_actual, u_control, sampleTime)
    
    current_eta = attitudeEuler(current_eta, nu, sampleTime)
    
    heading=current_eta[5]*180/math.pi %360
    
    speed=math.sqrt(nu[0]**2+nu[1]**2)
    # print('--------------------')

    # output=[current_eta[0],current_eta[1],u_actual[0],u_actual[1],heading,nu]
    output = [current_eta,nu,u_actual]
    
    # time.sleep(sampleTime) # simulasyon süresi ayarlama
    
    return output 

def control_allocation(u_avr,u_diff):
    max_frw_rpm = 104
    max_rvs_rpm = -104
    global u_control

    if u_avr>max_frw_rpm:
        u_avr=max_frw_rpm
    elif u_avr<max_rvs_rpm:
        u_avr=max_rvs_rpm
    # print('u_avr---',u_avr)
    n1=u_avr-u_diff
    n2=u_avr+u_diff

    if n1>max_frw_rpm:
        n1=max_frw_rpm
    if n1<max_rvs_rpm:
        n1=max_rvs_rpm
    if n2>max_frw_rpm:
        n2=max_frw_rpm   
    if n2<max_rvs_rpm:
        n2=max_rvs_rpm

    u_control=[n1,n2]
    
    return u_control

def heading_control(setpoint):
    global current_eta
    global heading
    global pre_error_heading 
    global saturation_limit_speed_max
    global saturation_limit_speed_min
    global prev_heading
    
    feed_back = heading
    error = setpoint-feed_back
    # if abs(error)<0.25:
    #     error=0
    print('error:',error)
    
    while error<-180:
        error = error +360
    while error > 180:
        error = error -360

    
    dt  = 0.02 
    ## PID  
    ## Proportion
    Up = -7*error

    ## Derivative
    # Ud = 10*(error-pre_error_heading)
    Ud = 0*(current_eta[-1]-prev_heading) # Aracın açısal konum farkından deribvative etkisi sisteme eklenir. Derivative kick etkisini azaltmak için.
    

    ## Integral 
    ki=125
    Ui =  ki* (error+pre_error_speed)*dt

    if Ui<saturation_limit_speed_max:
            Ui = Ui+  ki* (error+pre_error_speed)*dt/2 
    else :
        pass           

    ## control signal 
    u_diff = Up+Ud+Ui
    pre_error = error
    heading=current_eta[5]*180/math.pi %360
    print('pre error:',pre_error)
    pre_error_heading=error
    return u_diff
    # return 100

def speed_control(set_point):
    global nu
    global pre_error_speed
    global U_i
    global saturation_limit_speed_max
    global saturation_limit_speed_min
    global U_sat_old
    global prev_x
    global prev_y
    
    vehicle_velocity=math.sqrt(nu[0]**2+nu[1]**2)
    error=set_point-vehicle_velocity
    if set_point==0:
            set_point=0.0001
    else:
            pass
    dt = 0.02
    normal = 1
    if abs(error)<0.15:
        error=0
    ## FPID
    # Feedforward 
    U_f = (3.684*set_point**3-23.44*set_point**2+67.35*set_point+12.3)
    ## Proportion 
    U_p = 20*error  # 3 default nolmal0 =3  nolmal1=
    """Propetionu 12 yaptığımda control sinyali çıktısı değerinde aşım oluşmakta fakat hızlı bir oturma gerçekleşmektedir."""
    print('U_proportion',U_p)
    ## derivative 
    # U_d = 5*(error-pre_error_speed)   # 1 default 
    U_d = 0 * math.sqrt((prev_x-current_eta[0])**2+(prev_y-current_eta[1])**2) ## Konum farkları alınarak türev oluşturulur  
                                                                               ## Derivative kick etkisini en aza indirmek için.
    ## integral 
    U_i = 100* (error+pre_error_speed)*dt/2  # 0.25 default nolmal0 =0.3  nolmal1=
    """İntegrator değerini 1000 yaptığımızda kontrolcü sinyalınde aşım olmakta aynı zamanda hızlı bir şekilde ref. hıza oturmaktadır."""
    if normal ==1:

        if U_i<saturation_limit_speed_max:
            U_i = U_i+  0.25* (error+pre_error_speed)*dt/2 

        else :
            pass
        ## control signal
            
        u_avg = U_f+U_p+U_d+U_i
        pre_error_speed=error
        
        return u_avg
        

    elif normal==0:
        U_sat =U_sat_old + U_i
        U_sat_old=U_sat
        ## control signal 
        u_avg = U_f+U_p+U_d+U_i+U_sat
 
        ## Saturation
        U_filt=max(min(u_avg, saturation_limit_speed_max), saturation_limit_speed_min)
        if u_avg!=U_filt:
            U_sat_old = U_filt-U_p-U_f

        pre_error_speed=error


        return U_filt

def filtred_referans(pre_speed_ref,speed_ref,filtred_signal,delta_rate=0.03):

    if pre_speed_ref == filtred_signal:
        print('1')
        pre_speed_ref = filtred_signal
    else:
        if filtred_signal<speed_ref:
            print('2')
            filtred_signal = filtred_signal+delta_rate
            if filtred_signal>speed_ref:
                print('3')
                filtred_signal=speed_ref
        elif filtred_signal>speed_ref:
            print('444444444444444444444444444444444444444')
            filtred_signal = filtred_signal - delta_rate
            
            if filtred_signal<speed_ref:
                filtred_signal=filtred_signal

        else:
            filtred_signal=speed_ref

    

    return filtred_signal

def filtred_heading_referans(pre_speed_ref,speed_ref,filtred_signal):

    if speed_ref-filtred_signal<0:
        if (speed_ref-filtred_signal)%360<(filtred_signal-speed_ref):

            filtred_signal = filtred_signal+0.36
            if filtred_signal>0:

                if (-filtred_signal)%360>speed_ref:
                    filtred_signal=speed_ref
        elif  (speed_ref-filtred_signal)%360>(filtred_signal-speed_ref):

            filtred_signal = filtred_signal - 0.36
            if filtred_signal<0:
                if (filtred_signal%360)>speed_ref:
                    filtred_signal=speed_ref
    elif speed_ref-filtred_signal>0:
        if (speed_ref-filtred_signal)<(filtred_signal-speed_ref)%360:
            filtred_signal = filtred_signal+0.36
            if (filtred_signal)>0:
                if filtred_signal>speed_ref:
                    filtred_signal=speed_ref
        elif  (speed_ref-filtred_signal)>(filtred_signal-speed_ref)%360:

            filtred_signal = filtred_signal - 0.36
            if filtred_signal<0:

                if (filtred_signal%360)<speed_ref:
                        filtred_signal=speed_ref
        
        
    return filtred_signal


import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

if __name__ == "__main__": 
    start(0.02)
    # Wpx = [10,20,30,40]
    # Wpy = [10,20,10,20]
    Wpx = [10,20,30,20]
    Wpy = [10,20,30,40]
    # Wpx = [10,30,50,70]
    # Wpy = [10,20,10,40]
    ### non-index method
    # Wpx = [10,20,30,40]
    # Wpy = [10,10,10,10] 
    

    e,a,b,c,d = generate_curve(Wpx,Wpy,0.01)
    # Wpx = [10,20,20]
    # Wpy = [10,10,30]
    

    los2 = LineofSight()
    los  = LineofSightfirst()
    index = True
    if index==True:
        print('LOS 2 Method' )
        method=los2
    elif index==False:
        print('LOS Method')
        method=los
    R_cal=R_Calculator()
  
    x_list=[]
    y_list=[]
    d_list=[]
    init_heading=0.0
    r_list=[]

    prev_x=0
    prev_y=0
    
    prev_heading=0
    ds = 0.01 #[m]
    coeff = method.path_generate(Wpx,Wpy)

    v_list=[]
    head_list=[]
    timeOtter=[]
    ref_speed = []
    u_control_sancak = []
    u_control_iskele = []
    pervane_sancak = []
    pervane_iskele = []
    eta_list=[]
    Time=0
    chi_d_list = []
    actualheading_list = []
    x_closest_list=[]
    y_closest_list=[]
    x_los_list=[]
    y_los_list=[]
    x_closest=[]
    y_closest=[]
    ramp_signal=0
    Iramp_signal=0
    Sramp_signal=0
    Spp_u_control=0
    Ipp_u_control=0
    iramp_signal=0
    sramp_signal=0
    prev_u_control = [0,0]
    pp_u_control = 0#[0,0]
    prev_speed_ref=0.0
    prev_heading_ref=0.0
    filtred_signal=0.1
    heading_filtred_signal=0.1
    for i in range(4000):
        if index==True:
            print('LOS 2 Method' )
      
        elif index==False:
            print('LOS Method')
        
        Time+=0.02
        timeOtter.append(Time)
        time1=time.time()
        output=method.execute(nu,current_eta,Wpx,Wpy)
        refU=output['U_desired']
        refChi=math.degrees(output['chi_d'])%360
        filtred_signal= filtred_referans(prev_speed_ref,refU,filtred_signal)
        heading_filtred_signal= filtred_heading_referans(prev_heading_ref,refChi,heading_filtred_signal)

        U_diff = heading_control(heading_filtred_signal) # controlller
        path_curv_radius=R_cal.R_cal(current_eta)
        # U_desired = speed_generate(path_curv_radius,U_diff) # adaptive curve 
       
        u_avg = speed_control(filtred_signal)  # linear speed control 

        
        u_control=control_allocation(u_avg,U_diff) # adjust req. thrust for each thruster.
    
        output=function(u_control) # runnig dynamic model of usv (otter)


        v_list.append(np.sqrt(nu[0]**2+nu[1]**2))  # for plotting
        head_list.append(current_eta[5]*180/math.pi %360) # for plotting
        x_list.append(current_eta[0]) # for plotting 
        y_list.append(current_eta[1]) # for plotting 
        eta_list.append(current_eta)
        ref_speed.append(refU)
        u_control_sancak.append(u_control[0])
        u_control_iskele.append(u_control[1])
        x_closest_list.append(method.x_closest)
        y_closest_list.append(method.y_closest)
        x_los_list.append(method.x_los)
        y_los_list.append(method.y_los)
        pervane_iskele.append(output[2][1])
        pervane_sancak.append(output[2][0])
        chi_d_list.append(refChi)
        actualheading_list.append(heading)

        
        
        radius=R_cal.R_cal(current_eta)

        # method.los_simulation(current_eta,Wpx,Wpy,u_control)


        if method.x_los==method.x_closest and method.y_los == method.y_closest:

            break 
        
        

    # ani=los2.animate(x_list,y_list,head_list,Wpx,Wpy,x_closest_list,y_closest_list,x_los_list,y_los_list,u_control,i)
    # figure, ax = plt.subplots()
    # animation = FuncAnimation(figure,
    #                       func = ani,
    #                       frames = np.arange(0, 10, 0.1), 
    #                       interval = 10)
    # plt.show()
        
    print(u_control_sancak)        
    print('yyyyyyyyyy')  
    plt.figure()
    plt.plot(x_list,y_list)
    plt.plot(a,b)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.figure()
    plt.plot(timeOtter,v_list)
    plt.plot(timeOtter,ref_speed)
    plt.xlabel('Zaman')
    plt.ylabel('HIZ')
    plt.figure()
    plt.plot(timeOtter,u_control_sancak)
    plt.plot(timeOtter,pervane_sancak)
    plt.title('Sancak')
    plt.figure()
    plt.plot(timeOtter,u_control_iskele)
    plt.plot(timeOtter,pervane_iskele)
    plt.title('İskele')

    plt.figure()
    plt.plot(timeOtter,chi_d_list)
    plt.plot(timeOtter,actualheading_list)
    plt.title('Heading')

    plt.xlabel('Zaman')
    plt.ylabel('Kontrol Sinyali')
    plt.show()    

