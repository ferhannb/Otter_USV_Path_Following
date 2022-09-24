import numpy as np
from gnc import attitudeEuler
import math
import matplotlib.pyplot as plt
from otterPID import start,function



class TargetVehicle():

    def __init__(self,x=35,y=15,theta=90,v=3,omega=0):
        self.x = x
        self.y = y
        self.theta = theta
        self.v = v
        self.omega = omega
        self.predict_x = 0
        self.predict_y = 0
        self.x_list=[]
        self.y_list=[]
        self.theta_list=[]
        

    def KinematicTarget(self,dt=0.02):

        A = np.array([[math.cos(math.radians(self.theta)),0],[math.sin(math.radians(self.theta)),0],[0,1]])
        B = np.array([[self.v],[self.omega]])
        x_dot = A@B
        self.x_next = self.x + dt*x_dot[0]
        self.y_next = self.y + dt*x_dot[1]
        self.theta = self.theta + dt*x_dot[2]
        self.course =math.degrees(math.atan2(self.y_next-self.y,self.x_next-self.x))
        self.course = [(self.course + 360) % 360]
        self.x = self.x_next
        self.y = self.y_next

        self.x=self.x[0]
        self.y=self.y[0]
        self.course=self.course[0]
        U_speed = math.sqrt(x_dot[0]**2+x_dot[1]*2)

        return self.x,self.y,self.course,U_speed

    def prediction_movement(self,x,y,course,dt=0.02):


        self.x_list.append(x)
        self.y_list.append(y)
        self.theta_list.append(course)
   

        if len(self.x_list)==101:
            self.x_list.remove(self.x_list[0])
            self.y_list.remove(self.y_list[0])
            self.theta_list.remove(self.theta_list[0])
            

        x_dot_list = np.diff(self.x_list)
        len_x = len(x_dot_list)
        if len(x_dot_list)==0:
            len_x=1
        self.x_dot = (sum(x_dot_list)/(len_x))/dt

        y_dot_list = np.diff(self.y_list)
        len_y = len(y_dot_list)
        if len(y_dot_list)==0:
            len_y = 1
        self.y_dot = (sum(y_dot_list)/(len_y))/dt


        theta_dot_list = np.diff(self.theta_list)
        theta_len = len(theta_dot_list)
        if len(theta_dot_list)==0:
            theta_len = 1
        self.theta_dot = (sum(theta_dot_list)/(theta_len))/dt
        print('theta dot',self.theta_dot)


        predict_x_list=[]
        predict_y_list=[]
        predict_theta_list=[]



        U = math.sqrt(self.x_dot**2+self.y_dot**2)
        x_ =x
        y_ =y
        theta_dot=self.theta_dot

        theta_ =course

        for _ in range(1000):
            
            A = np.array([[math.cos(math.radians(theta_)),0],[math.sin(math.radians(theta_)),0],[0,1]])
            B = np.array([[U],[theta_dot]])
            x_dot = A@B

            x_next = x_ + dt*x_dot[0]
            y_next = y_ + dt*x_dot[1]
            theta_next= theta_ + dt*x_dot[2]
            theta_ =math.degrees(math.atan2(y_next-y_,x_next-x_))
            theta_ = [(theta_ + 360) % 360]

            x_ = x_next
            y_ = y_next
            theta_ = theta_next

            predict_x_list.append(x_next )
            predict_y_list.append(y_next)
            predict_theta_list.append(theta_)

        return predict_x_list,predict_y_list,predict_theta_list,theta_dot



    def simulation(self,predict_x_list,predict_y_list,x_o,y_o):

        plt.clf()
        plt.xlim([0,120])
        plt.ylim([0,120])
        plt.gcf().gca().add_artist(plt.Circle((self.x,self.y),2,fill=False))
        plt.plot(predict_x_list,predict_y_list,'r.-')
        plt.plot(x_o,y_o,'m.-')
        plt.show(block=False)
        plt.pause(0.02)


if __name__ == "__main__":
    targetVehicle = TargetVehicle() 
    otterVehicle = TargetVehicle()
    x=[]
    y=[]
    x_list=[]
    y_list=[]
    t_list=[]
    start(0.02)
    for i in range(1000):
        

        output = function([100,100])
        o_x,o_y,o_c = otterVehicle.prediction_movement(output[0][0],output[0][1],(math.degrees(output[0][5])%360))
        x_,y_,co,U_speed=targetVehicle.KinematicTarget()
        x_p,y_p,c_p=targetVehicle.prediction_movement(x_,y_,co)
        x_dist=[abs(a_x-b_x) for a_x, b_x in zip(o_x,x_p)]
        minvalue=min(x_dist)
        minindex=x_dist.index(minvalue)
        
        y_dist=[abs(a_y-b_y) for a_y, b_y in zip(o_y,y_p)]
        min_x =min(x_dist)
        min_y =min(y_dist)
        mindistV2V=math.sqrt(min_x**2+min_y**2)
        
        targetVehicle.simulation(x_p,y_p,o_x,o_y)


        