
import numpy as np
import math
from controllers.Controller import MPCController
from OtterDynamicClass import Otter
import matplotlib.pyplot as plt 
from kino_rrt.rrt_2D.kino_rrt2 import Kino_rrt
# from controllers.USVController import USVController
import random
import time 
 
### Path Smoothing 
def moving_average_filter(arr,window):
  
    window_size = window
    i = 0
    moving_averages = []
    while i < len(arr) - window_size + 1:
        
        # Store elements from i to i+window_size
        # in list to get the current window
        window = arr[i : i + window_size]
        # Calculate the average of current window
        window_average = round(sum(window) / window_size, 2)
        # Store the average of current
        # window in moving average list
        moving_averages.append(window_average)
        # Shift window to right by one position
        i += 1
    return moving_averages


window = 15
Wpx = []
Wpy = []
Xposelist = []
Yposelist = []
x_start = (5, 5, 0, 0, 0)  # Starting node
x_goal = (26.5, 19.5, math.pi/4, 0, 0)  # Goal node
# x_goal = (30, 23, math.pi/4, 0, 0)
MPC_otter = Otter()
MPC_otter.initialize_otter()
MPC_otter.current_eta=np.array([5,5,0,0,0,math.radians(30)]) ## Initial position of Otter Vehicle 

## Kino-RRT
errt = Kino_rrt(x_start, x_goal, 1, 0.5, 1, 20000, 1.5, 0, 0.35, -0.35, math.pi/2, -math.pi/2, 240, -140, 95, -95)
start = time.time()
errt.planning()
stop = time.time()
path_usv = errt.path[::-1]

path_usv_list = [list(ele) for ele in path_usv]
Xposelist = [sublist[0] for sublist in path_usv_list]
Yposelist = [sublistY[1] for sublistY in path_usv_list]
## Path Smoothing process
remain = len(Xposelist) // window
FiltXposelist = moving_average_filter(Xposelist,window)
FiltYposelist = moving_average_filter(Yposelist,window)
FiltXposelistL = FiltXposelist+Xposelist[-5:]
FiltYposelistL = FiltYposelist+Yposelist[-5:]
x_points = FiltXposelistL  ## Smoothed path X values
y_points = FiltYposelistL  ## Smoothed path Y values

## Determining points on genereted path
path_len_x = int(len(x_points))
path_len_y = int(len(y_points))
step_number = 200   ## Number of sub target point on genereted path
step_size_x = int(path_len_x/step_number)
step_size_y = int(path_len_y/step_number)
iter=0



for i in range(len(path_usv)):
    Wpx.append(path_usv[i][0])
    Wpy.append(path_usv[i][1])
    



init_states = [0,0,0,0,0,0]

HORIZON_STEP_COUNT = 30#100
MPC_DT = 0.1 ## sapling time of MPC cost function
MPC_REFRESH_TIME = 0.1# how many horizon time we will wait for recalculating mpc

Q30 = [100,0,760,760,0] #[1000,0,260,260,0]
R30 = [100,100]
LOS_MPC = MPCController(N=HORIZON_STEP_COUNT, dt=MPC_DT,Q=Q30,R=R30)

SIM_DT = 0.02
SIM_TIME = 150*13
time_ = 0


time_list = []
#### MPC lists
iskeleL=[]
sancakL=[]
time_l = []
x_list=[]
y_list=[]
port_rpm_list=[]
starboard_rpm_list=[]
timer=0
rand_speed_dist = random.uniform(-0.05,0.05)
inlet_x_path=[]
inlet_y_path=[]
x_target_list = []
y_target_list = []
mpc_x = []
mpc_y = []
c = 0
for i in range(SIM_TIME):

    timer+=0.02

    if iter<=(step_number-1):      
        print('dsd')
        target_x = x_points[iter*step_size_x]
        target_y = y_points[iter*step_size_y]
        target_x_next = x_points[(iter+1)*step_size_x-1]
        target_y_next = y_points[(iter+1)*step_size_y-1]
    elif iter==step_number:
        target_x = x_points[-2]
        target_y = y_points[-2]
        target_x_next = x_points[-1]
        target_y_next = y_points[-1]

    x_target_list.append(target_x)
    y_target_list.append(target_y)
    dx = target_x_next-target_x
    dy = target_y_next-target_y
    target_psi =  math.atan2(dy, dx)
    target_values=[1, 0, 0, target_x, target_y, target_psi]

    print("ITERATION NUMBER:",iter)
    print('path length', path_len_x)
    print("Target X {} - Position X {}".format(target_x,MPC_otter.current_eta[0]))
    print("Target Y {} - Position Y {}".format(target_y,MPC_otter.current_eta[1]))
    print("Target Psi {} - Actual Psi {} ".format(math.degrees(target_psi),math.degrees(MPC_otter.current_eta[-1])))
    print('TARGET PSİ',math.degrees(target_psi))
    
    if math.sqrt((target_x-MPC_otter.current_eta[0])**2+(target_y-MPC_otter.current_eta[1])**2)<2:   
        print('close') 
        iter += 1
        
    time_list.append(timer)



                                ### DISTURBANCE ####

    # if i<=40/SIM_DT:
    #     beta=0
    #     vc=0.1
    #     MPC_otter.current_disturbance(vc+rand_speed_dist,beta)
    #     mpcpid_otter.current_disturbance(vc+rand_speed_dist,beta)
    #     pidOtter.beta_c =beta
    #     pidOtter.V_c = vc+rand_speed_dist
    #     inletOtter.beta_c = beta
    #     inletOtter.V_c = vc+rand_speed_dist        
    # elif i<=80/SIM_DT:
    #     beta=56
    #     vc = 0.3
    #     MPC_otter.current_disturbance(vc+rand_speed_dist,beta)
    #     # mpcpid_otter.current_disturbance(0.1+rand_speed_dist,beta)
    #     pidOtter.beta_c =beta
    #     pidOtter.V_c = vc+rand_speed_dist
    #     inletOtter.beta_c = beta
    #     inletOtter.V_c = vc+rand_speed_dist     
    # else:
    #     beta=124
    #     vc=0.2
    #     MPC_otter.current_disturbance(vc+rand_speed_dist,beta)
    #     mpcpid_otter.current_disturbance(vc+rand_speed_dist,beta)
    #     pidOtter.beta_c =beta
    #     pidOtter.V_c = vc+rand_speed_dist
    #     inletOtter.beta_c = beta
    #     inletOtter.V_c = vc+rand_speed_dist


  

    if time_ == SIM_DT:
        
        dt_cnt = 0
        OneMPC_output=LOS_MPC.execute_MPC(cnt = i, init_states=[MPC_otter.nu[0],MPC_otter.nu[5],MPC_otter.current_eta[0],MPC_otter.current_eta[1],MPC_otter.current_eta[5]],target_states=[0.7,0,target_x,target_y,target_psi])
    elif  i % int(MPC_REFRESH_TIME / SIM_DT) == 0: # MPC'nin execute time 
        OneMPC_output=LOS_MPC.execute_MPC(cnt = i, init_states=[MPC_otter.nu[0],MPC_otter.nu[5],MPC_otter.current_eta[0],MPC_otter.current_eta[1],MPC_otter.current_eta[5]],target_states=[0.7,0,target_x,target_y,target_psi])
        dt_cnt = 0
    if i % int(MPC_DT / SIM_DT) == 0:
        iskele_rpm = MPC_otter.thruster(OneMPC_output['input_prediction'][0][dt_cnt])
        sancak_rpm = MPC_otter.thruster(OneMPC_output['input_prediction'][1][dt_cnt])
        MPC_otter.u_control = [iskele_rpm,sancak_rpm]
        dt_cnt += 1

    output_func=MPC_otter.function(sampleTime=SIM_DT)
    c+=1
   
    x_init=output_func['current_eta'][0]
    y_init=output_func['current_eta'][1]
    yaw_init=output_func['current_eta'][5]

    mpc_x.append(MPC_otter.current_eta[0])
    mpc_y.append(MPC_otter.current_eta[1])
    port_rpm_list.append(MPC_otter.u_control[1])
    starboard_rpm_list.append(MPC_otter.u_control[0])
        


    # if np.sqrt((FiltXposelistL[-1]-x_init)**2+(FiltYposelistL[-1]-y_init)**2)<1:
    #     break

## Scenario-1 

from matplotlib.patches import Rectangle
fig, ax = plt.subplots()
scenario = 2
if scenario ==1:
    ax.add_patch(
                    Rectangle(
                    (14, 12), 4, 2,
                    edgecolor='black',
                    facecolor='black',
                    fill=True
                    )
                )

    ax.add_patch(Rectangle(
                    (18, 22), 4, 3,
                    edgecolor='black',
                    facecolor='black',
                    fill=True
                    )
                )

    ax.add_patch(Rectangle(
                    (26, 7), 2, 4,
                    edgecolor='black',
                    facecolor='black',
                    fill=True
                    )
                )
    ax.add_patch(Rectangle(
                    (32, 14), 4, 2,
                    edgecolor='black',
                    facecolor='black',
                    fill=True
                    )
                )
## Scenario-2
if scenario ==2:
    ax.add_patch(
                    Rectangle(
                    (22, 15), 2, 6,
                    edgecolor='black',
                    facecolor='black',
                    fill=True
                    )
                )

    ax.add_patch(Rectangle(
                    (24, 15), 6, 2,
                    edgecolor='black',
                    facecolor='black',
                    fill=True
                    )
                )

    ax.add_patch(Rectangle(
                    (29, 15), 2, 6,
                    edgecolor='black',
                    facecolor='black',
                    fill=True
                    )
                )
    ax.add_patch(Rectangle(
                    (4, 15), 5, 3,
                    edgecolor='black',
                    facecolor='black',
                    fill=True
                    )
                )

ax.add_patch(
                Rectangle(
                (12, 5), 1, 3,
                edgecolor='black',
                facecolor='black',
                fill=True
                )
            )



ax.plot(FiltXposelistL,FiltYposelistL,linewidth='2',label='Generated smoothed Path')
ax.plot(mpc_x,mpc_y,linewidth='2',label = 'MPC Path Following Performance')
ax.plot(FiltXposelistL[0],FiltYposelistL[0],"g o",label='Start Point')
ax.plot(FiltXposelistL[-1],FiltYposelistL[-1],"r o",label='Goal Point')

plt.text(FiltXposelistL[0],FiltYposelistL[0],'Start Point',fontsize=10,bbox=dict(facecolor='red',alpha=0.7))
plt.text(FiltXposelistL[-1],FiltYposelistL[-1],'Goal Point',fontsize=10, bbox=dict(facecolor='red',alpha=0.7))
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.title('Path Following Controller Performance')
plt.legend()
print('X',FiltXposelist[0])
print('Y',FiltYposelistL[0])
plt.savefig('Seneryo1.png',dpi=300)
plt.show()

# plt.figure()
# plt.plot(time_list,sancak_rpm_list, label='Port RPS')
# plt.plot(time_list,iskele_rpm_list,label='Starboard RPS')

# plt.title('MPC Input Values')
# plt.legend()
# plt.xlabel('TIME (s)')
# plt.ylabel('RPS')
# plt.savefig('DistYeniSc2MPCIputgraphSc1.png',dpi=300)
# plt.show()


plt.figure()
plt.plot(time_list,port_rpm_list,label='Port RPS')
plt.plot(time_list,starboard_rpm_list, label='Starboard RPS')
plt.title('MPC Input Values')
plt.legend()
plt.xlabel('TIME (s)')
plt.ylabel('RPS')
plt.savefig('DistYeniSc2FPIDIputgraph.png',dpi=300)
plt.show()
