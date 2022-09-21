#!/usr/bin/env python
# -*- coding: utf-8 -*-


from USVController import USVController
from OtterDynamicClass import Otter
from path_track_2 import LineofSight
from path_track_3 import LineofSightfirst
import matplotlib.pyplot as plt 

pid = USVController()
vehicle = Otter()
path_generation = LineofSight()
path_generation_diff = LineofSightfirst()




vehicle.initialize_otter()
x_list = []
y_list = []
speed_list = []
heading_list =  []
time_list = []
time_=0
Kp_heading = -10
Ki_heading = -4
Kd_heading = 0
Kp_speed = 1
Ki_speed = 1
Kd_speed = 0

for _ in range(4000):
    time_+=0.02
    ref_heading = 60 
    ref_speed = 3
    U_diff=pid.Heading_controller(ref_heading,vehicle.current_eta,Kp_heading,Ki_heading,Kd_heading)
    # u_avg = pid.Speed_controller(ref_speed)
    vehicle.u_control = pid.control_allocation(0,U_diff)
    output = vehicle.function()
    print(output['heading'])
    speed_list.append(output['speed'])
    time_list.append(time_ )
    x_list.append(vehicle.current_eta[0])
    y_list.append(vehicle.current_eta[1])
    heading_list.append(vehicle.current_eta[-1])
plt.plot(time_list,heading_list)
plt.show()
plt.plot(x_list,y_list)
plt.show()


