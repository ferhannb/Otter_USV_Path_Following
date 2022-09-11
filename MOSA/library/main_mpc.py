
import numpy as np
import math
import time
import pandas as pd
from gnc import attitudeEuler
from path_track_2 import LineofSight
from R_calculation import R_Calculator
from mpc_colav import MPCUSV_colav
from mpccontroller import MPCUSV
import casadi as ca 
import time
from otterMPC import Otter
from dynamic_obs import TargetVehicle
from colreg_define import COLREG_detect,colreg_situation


### akıntı hızı ve yönü 
V_current = 0. # akıntı hızı
beta_current =0 # akıntı yönü

def start(steptime):
    # vehicle=Otter()
    global current_eta
    global nu
    global u_actual
    global u_control
    global sampleTime
    global speed
    global heading

    
    nu = np.array([3, 0, 0, 0, 0, 0], float)
    u_actual = np.array([0, 0], float)
    u_control = np.array([0, 0], float)           
    current_eta = np.array([10, 10, 0, 0, 0, math.radians(60)], float)
    heading=current_eta[5]*180/math.pi %360
    speed = math.sqrt(nu[0]**2+nu[1]**2)
    sampleTime=steptime



def function(u_control):

    vehicle=Otter(V_current = V_current, 
                beta_current = beta_current,)

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




#####################################################
#################### SIMULATION #####################
#####################################################


start(0.02)
# Wpx = [0,100]
# Wpy = [30,30]
Wpx = [10,20,30]
Wpy = [10,20,30]
# Wpx = [10,45,60,80]
# Wpy = [20,30,70,30]
los2=LineofSight()
R_cal=R_Calculator()

# obstacle = TargetVehicle(40,10,90,2,0) # SAFE   #(40,10,90,2,0) # Give away    #(40,10,90,2,-6)  çember
obstacle = TargetVehicle(-100,-100,90,0.)
OtterPredict = TargetVehicle(0,0,90,0)

def set_obstacle():
    """set_obstacle fonksiyonu dinamik engel oluştrumak için kullanilir. (X,Y,heading,velocity) parametrelerini alabilir."""

    diam_usv = 2
    diam_obs = 3
    x_obs,y_obs,course,obs_speed=obstacle.KinematicTarget()   
    x_obsArr,y_obsArr,courseArr,TS_theta_dot=obstacle.prediction_movement(x_obs,y_obs,course)

    dist = math.sqrt((current_eta[0]-x_obs)**2+(current_eta[1]-y_obs)**2)-(diam_usv/2+diam_obs/2+0.5) 
    safe_dist = 25

    if dist>safe_dist:
        param=True

    else:
        param=False
    return param,x_obs,y_obs,safe_dist,course,x_obsArr,y_obsArr,courseArr,obs_speed,TS_theta_dot


def colreg_execute(param,init_states,los_output,colreg):


    if param:
        # print('LOS')
        mpc_output= mpc_los.execute_MPC(init_states,[los_output['U_desired'],0,0,los_output['x_los'],los_output['y_los'],los_output['chi_d']])
    else:
        print('==OBSTACLE')
        if colreg['COLREG']=='Head-on':
            # print('HEAD_ON')
            sc='HEAD-ON'
            mpc_output= mpc.execute_MPC(init_states,[los_output['U_desired'],0,0,los_output['Wpx'],los_output['Wpy'],math.atan2((los_output['Wpy']-current_eta[1]),(los_output['Wpx']-current_eta[0]))])
            if 180<=colreg['RelTarget'] and colreg['RelTarget']<=270:
                mpc_output= mpc_los.execute_MPC(init_states,[0,0,0,los_output['x_los'],los_output['y_los'],los_output['chi_d']])

        elif colreg['COLREG']=='Give-Away':
            # print('Give-Away')
            mpc_output= mpc.execute_MPC(init_states,[los_output['U_desired'],0,0,x_obs+col_dist*math.cos(math.radians(course+180))
                        ,y_obs+col_dist*math.sin(math.radians(course+180)),math.atan2((y_obs-current_eta[1]),(x_obs+10-current_eta[0]))])
            # mpc_output= mpc.execute_MPC(init_states,[0,0,0,x_obs+10,y_obs,math.atan2((y_obs-current_eta[1]),(x_obs+10-current_eta[0]))])
            # mpc_output= mpc.execute_MPC(init_states,[0,0,0,los_output['Wpx'],los_output['Wpy'],math.atan2((los_output['Wpy']-current_eta[1]),(los_output['Wpx']-current_eta[0]))])
        elif colreg['COLREG']=='Safe':
            # print('11Safe')
            mpc_output= mpc_los.execute_MPC(init_states,[los_output['U_desired'],0,0,los_output['x_los'],los_output['y_los'],los_output['chi_d']])
    
    return colreg,mpc_output



    
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
init_states = [0,0,0,30,30,0]


for i in range(8000):
    
    
    los_output=los2.execute(nu,current_eta,Wpx,Wpy)
    
    param,x_obs,y_obs,safe_dist,course,x_obsArr,y_obsArr,courseArr,TS_U_speed,theta_dot=set_obstacle()

    TS_Arr=[x_obsArr,y_obsArr,courseArr]

    mpc_los = MPCUSV()
    mpc = MPCUSV_colav(course,x_obs,y_obs)

    
    ship1 = pd.DataFrame({'MMSI':[1], 'LAT':current_eta[0], 'LON':current_eta[1], 'SOG':math.sqrt(nu[0]**2+nu[1]**2), 'COG':math.degrees(current_eta[5])})
    ship2 = pd.DataFrame ({'MMSI':[1], 'LAT':x_obs, 'LON':y_obs, 'SOG':obstacle.v, 'COG':course})
    OS_xArr,OS_yArr,OS_courseArr,OS_theta_dot=OtterPredict.prediction_movement(current_eta[0],current_eta[1],current_eta[5])

    OS_Arr =[OS_xArr,OS_yArr,OS_courseArr]
    colreg = COLREG_detect(ship1,ship2) 
    # print('COLREG DETECT',colreg)
    col_dist  =5
    OS_U_speed = math.sqrt(nu[0]**2+nu[1]**2)
    colreg = colreg_situation(colreg,OS_U_speed,current_eta,OS_xArr,OS_yArr,OS_courseArr,TS_U_speed,x_obs,y_obs,course,x_obsArr,y_obsArr,courseArr,theta_dot)
    # print('COLREG SİTUATION',colreg)
    colreg,mpc_output = colreg_execute(param,init_states,los_output,colreg)
    
  


    u_d = np.array(mpc_output['control_signal'].full())
    u_control = [float(u_d[0][0]),float(u_d[1][0])]


    v_list.append(np.sqrt(nu[0]**2+nu[1]**2))
    head_list.append(current_eta[5]*180/math.pi %360)
    x_list.append(current_eta[0])
    y_list.append(current_eta[1])

    output=function(u_control)
    x_init=output[0][0] 
    y_init=output[0][1]
    yaw_init=output[0][5]

    init_states = ca.vertcat(0,0,0,x_init,y_init,yaw_init)
    radius=R_cal.R_cal(current_eta)
    los2.los_simulation(current_eta,Wpx,Wpy,u_control,x_obs,y_obs,[col_dist*math.cos(math.radians(course+180)+col_dist*math.sin(math.radians(course+180))),colreg,],OS_Arr,TS_Arr)