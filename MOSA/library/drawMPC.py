
import numpy as np
import math
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches



class Draw_MPC_point_stabilization(object):

    def __init__(self, state_history: np.array,usv_states: list, init_state: np.array, target_state: np.array,N=30,export_fig=False):

        self.usv_states = usv_states
        self.init_states = init_state
        self.target_states = target_state
        self.trajectory = state_history
        self.N=N
        self.diam = 2
        self.fig = plt.figure()
        self.ax = plt.axes(xlim = (-25,25),ylim=(-25,25))
        self.fig.set_size_inches(7,7)
        self.animation_init()

        self.ani = animation.FuncAnimation(self.fig, func=self.animation_loop, 
                                           init_func=self.animation_init, frames=len(self.usv_states),interval=100, repeat=False)

        plt.grid('--')
        plt.show()
        
        

    def animation_init(self):
 
        ## usv states  
        
        x_usv     = int(self.usv_states[0][3])
        y_usv     = int(self.usv_states[0][4])
        theta_usv = int(self.usv_states[0][5]) 
        # x_usv     = int(self.usv_states[0][0])
        # y_usv     = int(self.usv_states[0][1])
        # theta_usv = int(self.usv_states[0][2]) 
        ## init states
        print(self.init_states)
        x_init = int(self.init_states[3])
        y_init = int(self.init_states[4])
        theta_init = int(self.init_states[5])

        # x_init = int(self.init_states[0])
        # y_init = int(self.init_states[1])
        # theta_init = int(self.init_states[2])
        ## target states 
        # x_tar = int(self.target_states[0])
        # y_tar = int(self.target_states[1])
        # theta_tar = int(self.target_states[2])

        x_tar = int(self.target_states[3])
        y_tar = int(self.target_states[4])
        theta_tar = int(self.target_states[5])
        
        traj_x = self.trajectory[3,:,0]
        traj_y = self.trajectory[4,:,0]

        # traj_x = self.trajectory[0,:,0]
        # traj_y = self.trajectory[1,:,0]

        self.target_pos = plt.Circle([x_tar,y_tar], self.diam, color='b', fill=False)
        self.ax.add_artist(self.target_pos)
        self.target_arr = mpatches.Arrow(x_tar, y_tar,
                                         self.diam * np.cos(theta_tar),
                                         self.diam * np.sin(theta_tar), width=0.4)

        self.ax.add_patch(self.target_arr)
        self.usv = plt.Circle([x_init,y_init], self.diam, color='g', fill=False)
        self.ax.add_artist(self.usv)
        self.usv_arr = mpatches.Arrow(x_init, y_init,
                                        self.diam * np.cos(theta_init),
                                        self.diam * np.sin(theta_init), width=0.4, color='k')

        self.ax.add_patch(self.usv_arr)
        self.horizon,=self.ax.plot(traj_x,traj_y,'-*r')

        return self.target_pos, self.target_arr, self.usv, self.usv_arr



    def animation_loop(self,idx):

        traj_x = self.trajectory[3,:,idx]
        traj_y = self.trajectory[4,:,idx]
        
        position = self.usv_states[idx][3:5]
        orientation = self.usv_states[idx][5]

        # traj_x = self.trajectory[0,:,idx]
        # traj_y = self.trajectory[1,:,idx]
        
        # position = self.usv_states[idx][0:2]
        # orientation = self.usv_states[idx][2]

        self.usv.center = position
        self.usv_arr.remove()
        self.usv_arr = mpatches.Arrow(position[0], position[1], self.diam * np.cos(orientation),
                                        self.diam * np.sin(orientation), width=0.2, color='r')
        self.ax.add_patch(self.usv_arr)
        
        self.horizon.remove()
        self.horizon,=self.ax.plot(traj_x,traj_y,'-.r')
        
        return self.usv_arr, self.usv, 
    

if __name__ =='__main__':

    draw =Draw_MPC_point_stabilization([],[0,0,math.pi/8],[-2,-2,0],[4,4,math.pi/2])
    draw.animation_init()