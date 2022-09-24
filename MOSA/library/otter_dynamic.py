import numpy as np
import math
import time


#from OtterSimulator.control import PIDpolePlacement
#from OtterSimulator.gnc import Smtrx, Hmtrx, m2c, crossFlowDrag, sat, attitudeEuler
from gnc import Smtrx, Hmtrx, m2c, crossFlowDrag, sat, attitudeEuler
from path_track_2 import LineofSight
from R_calculation import R_Calculator
from otterradius import speed_generate
from mpccontroller import MPCUSV
import casadi as ca 
import time
import scipy as sp
from otter2 import Otter
def start(steptime):
    vehicle=Otter()
    global current_eta
    global nu
    global u_actual
    global u_control
    global sampleTime
    global speed
    global heading

    
    nu = np.array([0, 0, 0, 0, 0, 0], float)
    u_actual = np.array([0, 0], float)
    u_control = np.array([0, 0], float)           
    current_eta = np.array([0, 10, 0, 0, 0, math.radians(0)], float)
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
    [nu, u_actual]=vehicle.dynamics(current_eta, nu, u_actual, u_control, sampleTime)
    current_eta = attitudeEuler(current_eta, nu, sampleTime)
    heading=current_eta[5]*180/math.pi %360
    speed=math.sqrt(nu[0]**2+nu[1]**2)
    # print('--------------------')

    # output=[current_eta[0],current_eta[1],u_actual[0],u_actual[1],heading,nu]
    output = [current_eta,nu,u_actual]
    time.sleep(sampleTime) # simulasyon süresi ayarlama
    return output 

def control_allocation(u_avr,u_diff):
    max_frw_rpm = 100
    max_rvs_rpm = -100
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

    feed_back = heading
    error = setpoint-feed_back

    while error<-180:
        error = error +360
    while error > 180:
        error = error -360

    Up = -5*error
    u_diff = Up

    heading=current_eta[5]*180/math.pi %360
    

    return u_diff
    # return 100

def speed_control(set_point):

    global nu
    vehicle_velocity=math.sqrt(nu[0]**2+nu[1]**2)
    error=set_point-vehicle_velocity
    if set_point==0:
            set_point=0.0001
    else:
            pass
    U_f = (3.684*set_point**3-23.44*set_point**2+67.35*set_point+12.3)
    U_p = 10*error
    u_avg = U_f+U_p
    # print('u_avg',u_avg)
    return u_avg
    



if __name__ == "__main__": 

    start(0.02)
    Wpx = [0,30,50]
    Wpy = [10,50,25]
    los2=LineofSight()
    R_cal=R_Calculator()
    mpc = MPCUSV()
    x_list=[]
    y_list=[]
    d_list=[]
    init_heading=0.0
    r_list=[]

    prev_x=0
    prev_y=0
    
    prev_heading=0
    ds = 0.01 #[m]
    coeff = los2.path_generate(Wpx,Wpy)
    v_list=[]
    head_list=[]
    init_states = [0,0,0,0,10,0]
    for i in range(4000):
        

        error_angle, U_desired,y_e,chi_d,target_angle,x_los,y_los=los2.execute(current_eta,Wpx,Wpy)
        state_history,xx,t,init_states,target_states,home_pose,u = mpc.execute_MPC(init_states,[0,0,0,x_los,y_los,target_angle])

        u_d = np.array(u.full())
        u_control = [float(u_d[0][0]),float(u_d[1][0])]


        v_list.append(np.sqrt(nu[0]**2+nu[1]**2))
        head_list.append(current_eta[5]*180/math.pi %360)
        x_list.append(current_eta[0])
        y_list.append(current_eta[1])
        # u_control=[-70,100]
        output=function(u_control)
        x_init=output[0][0] 
        y_init=output[0][1]
        yaw_init=output[0][5]

        init_states = ca.vertcat(0,0,0,x_init,y_init,yaw_init)
        radius=R_cal.R_cal(current_eta)
        los2.los_simulation(current_eta,Wpx,Wpy,u_control)
   
