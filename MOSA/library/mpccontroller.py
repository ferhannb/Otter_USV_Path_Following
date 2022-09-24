

import casadi as ca 
import numpy as np
import math 
import time
from otterMPC import Otter
from drawMPC import Draw_MPC_point_stabilization



class MPCUSV():

    def __init__(self, N=70, dt = 0.15):

        self.N = N
        self.dt = dt

        self.eta = [0,0,0,0,0,0]
        
        ## states 

        u = ca.SX.sym('u')
        v = ca.SX.sym('v')
        r = ca.SX.sym('r')
        x = ca.SX.sym('x')
        y = ca.SX.sym('y')
        psi = ca.SX.sym('psi')   
        states=ca.vertcat(u,v,r,x,y,psi)
        self.n_states =states.numel()

        # inputs otter

        Tp = ca.SX.sym('Tp')
        Ts = ca.SX.sym('Ts')
        inputs = ca.vertcat(Tp,Ts)
        B = 1.08/2 # Beam  
        self.n_inputs = inputs.numel()
        # tau = ca.vertcat(Tp+Ts,0,(B*(Tp-Ts)))

        otter = Otter()
        self.D,self.C,self.M = otter.dynamic_eq()

      

        dynamics = ca.vertcat((self.M[0][0]+self.M[1][0]+self.M[2][0])*((Tp+Ts)+(self.D[0][0]+self.D[1][0]+self.D[2][0])*u+(self.C[0][0]+self.C[1][0]+self.C[2][0])*u),
                              (self.M[0][1]+self.M[1][1]+self.M[2][1])*((0)+(self.D[0][1]+self.D[1][1]+self.D[2][1])*v+(self.C[0][1]+self.C[1][1]+self.C[2][1])*v),
                              (self.M[0][2]+self.M[1][2]+self.M[2][2])*(B*(Tp-Ts)+(self.D[0][2]+self.D[1][2]+self.D[2][2])*r+(self.C[0][2]+self.C[1][2]+self.C[2][2])*r),
                              u*ca.cos(psi)-v*ca.sin(psi),
                              u*ca.sin(psi)+v*ca.cos(psi),
                              r)

        # dynamics = ca.vertcat ((0.01262154610383696)*((Tp+Ts)+(77.55443234836703)*u+(0)*u),
        #                        (0.006560008561141038)*(0+(-0.0)*v+(0)*v),
        #                        (0.00799364586576341)*(B*(Tp-Ts)+(605.4038837833798)*r+(0)*r),
        #                         u*ca.cos(psi)-v*ca.sin(psi),
        #                         u*ca.sin(psi)+v*ca.cos(psi),
        #                         r)
                                

              
        Q = ca.SX.zeros(self.n_states,self.n_states)
        Q[0,0] = 1#500
        Q[1,1] = 1#0.0
        Q[2,2] = 10000#100000
        Q[3,3] = 30#3500
        Q[4,4] = 30#3500
        Q[5,5] = 10000#1000

        QF = ca.SX.zeros(self.n_states,self.n_states)
        QF[0,0] = 1
        QF[1,1] = 1
        QF[2,2] = 100000
        QF[3,3] = 40000
        QF[4,4] = 40000
        QF[5,5] = 100000

        # Q[0,0] = 0.01#500
        # Q[1,1] = 0.01#0.0
        # Q[2,2] = 100#1000000#1
        # Q[3,3] = 100#350#250
        # Q[4,4] = 100#150#44#300
        # Q[5,5] = 1000#300000#3000#500
        
        
        R = ca.SX.zeros(self.n_inputs,self.n_inputs)
        R[0,0] = 0.025
        R[1,1] = 0.025




        self.f = ca.Function('f',[states,inputs],[dynamics])
        
        # Optimization variables for all states and inputs in control horizon


        X = ca.SX.sym('X',self.n_states,self.N+1)
        U = ca.SX.sym('U',self.n_inputs,self.N)
        P = ca.SX.sym('P',self.n_states+self.n_states)

        cost_function = 0
        
        g = X[:, 0] - P[:self.n_states] 
        # g = []
        # g.append([X[:, 0] - P[:self.n_states]])

        for k in range(0,self.N-1):

            st = X[:,k]
            input =U[:,k]
            cost_function = cost_function + (st-P[self.n_states:]).T@ Q @ (st-P[self.n_states:])+ input.T @ R @ input
            st_next = X[:,k+1]
            f_value = self.f(st,input)
            st_next_euler = st + (self.dt*f_value)
            g = ca.vertcat(g,st_next-st_next_euler)
            # g.append([st_next-st_next_euler])

        cost_function=cost_function+(X[:,self.N]-P[self.n_states:]).T@QF@(X[:,self.N]-P[self.n_states:])
        
        # contain all decision variable into one row vector

        # OPT_variables = ca.vertcat(ca.reshape(X,self.n_states*(self.N+1),1),ca.reshape(U,self.n_inputs*(self.N),1))

        OPT_variables = ca.vertcat( ca.reshape(X, -1, 1), ca.reshape(U, -1, 1))

        nlp_prob = {
                    'f': cost_function,
                    'x': OPT_variables,
                    'g': g,
                    'p': P
                    }


        opts = {
                'ipopt': {
                    'max_iter': 2000,
                    'print_level': 0,
                    'acceptable_tol': 1e-8,
                    'acceptable_obj_change_tol': 1e-6
                },
                'print_time': 0
            }
            


        self.solver = ca.nlpsol('solver','ipopt',nlp_prob,opts)


        lbx = ca.DM.zeros((self.n_states*(self.N+1) + self.n_inputs*self.N, 1))
        ubx = ca.DM.zeros((self.n_states*(self.N+1) + self.n_inputs*self.N, 1))





        lbx[0: self.n_states*(N+1): self.n_states] = -3         # u lower bound
        lbx[1: self.n_states*(N+1): self.n_states] = -0.5         # v lower bound
        lbx[2: self.n_states*(N+1): self.n_states] = -0.35     # r lower bound
        lbx[3: self.n_states*(N+1): self.n_states] = -ca.inf     # x lower bound
        lbx[4: self.n_states*(N+1): self.n_states] = -ca.inf     # y lower bound
        lbx[5: self.n_states*(N+1): self.n_states] = -ca.inf     # psi lower bound

        ubx[0: self.n_states*(N+1): self.n_states] = 3         # u lower bound
        ubx[1: self.n_states*(N+1): self.n_states] = 0.5         # v lower bound
        ubx[2: self.n_states*(N+1): self.n_states] = 0.35     # r lower bound
        ubx[3: self.n_states*(N+1): self.n_states] = ca.inf     # x lower bound
        ubx[4: self.n_states*(N+1): self.n_states] = ca.inf     # y lower bound
        ubx[5: self.n_states*(N+1): self.n_states] = ca.inf     # psi lower bound
 
        lbx[self.n_states*(self.N+1)::self.n_inputs] =   -100    # Tp lower bound
        ubx[self.n_states*(self.N+1)::self.n_inputs] =    100    # Tp upper bound
        lbx[self.n_states*(self.N+1)+1::self.n_inputs] = -100 # Ts lower bound
        ubx[self.n_states*(self.N+1)+1::self.n_inputs] =  100 # Tp upper bound
        
        
        
        # equality constraints 

        lbg = ca.DM.zeros(self.n_states*(self.N),1)
        ubg = ca.DM.zeros(self.n_states*(self.N),1)
       
        # lbg[0:self.n_states] = 0.0
        # ubg[0:self.n_states] = 0.0

        # self.args = {
        #         'lbg': ca.DM.zeros((self.n_states*(N+1), 1)),  # constraints lower bound
        #         'ubg': ca.DM.zeros((self.n_states*(N+1), 1)),  # constraints upper bound
        #         'lbx': lbx,
        #         'ubx': ubx
        #
        #              }
        self.args = {
                'lbg': lbg,  # constraints lower bound
                'ubg': ubg,  # constraints upper bound
                'lbx': lbx,
                'ubx': ubx
                    }

    def shift_timestep(self,dt, t0, state_init, u, f):

        f_value = f(state_init, u[:, 0])
        next_state = ca.DM.full(state_init + (dt * f_value))

        t0 = t0 + dt
        u0 = ca.horzcat(
            u[:, 1:],
            ca.reshape(u[:, -1], -1, 1))

        return t0, next_state, u0


    def DM2Arr(self,dm):
        return np.array(dm.full())
    

    def execute_MPC(self,init_states=[0,0,0,0,0,math.pi/2],target_states=[0,0,0,10,10,0]):

        init_states = np.array(init_states)
        home_pose = init_states
        
        # target_states = np.array([x_los,y_los,heading])
        t0=0

        t = ca.DM(t0)
        u0 = ca.DM.zeros((self.N,self.n_inputs))
        X0 = ca.repmat(init_states,1,self.N+1)
        
        # MPC information
        mpc_iteration = 0
        sim_time = 40

        state_history = self.DM2Arr(X0)
        input_history = self.DM2Arr(u0[:, 0])
        xx=[]
        times = np.array([[0]])
        
        u_c=[] 

        # while (ca.norm_2(init_states - target_states) > 1e-1) and (mpc_iteration * self.dt < sim_time):
        for i in range(1):
            t1 = time.time()
            self.args['p'] = ca.vertcat(
                init_states,    # current state
                target_states   # target state
            )
            # optimization variable current state
            self.args['x0'] = ca.vertcat(
                ca.reshape(X0, self.n_states*(self.N+1), 1),
                ca.reshape(u0, self.n_inputs*self.N, 1)
            )

            sol = self.solver(
                x0  = self.args['x0'],
                lbx = self.args['lbx'],
                ubx = self.args['ubx'],
                lbg = self.args['lbg'],
                ubg = self.args['ubg'],
                p   = self.args['p']
                             )

            u  = ca.reshape(sol['x'][self.n_states * (self.N + 1):], self.n_inputs, self.N)
            u_c.append(u[:,0])
            X0 = ca.reshape(sol['x'][: self.n_states * (self.N+1)], self.n_states, self.N+1)
            
            t0, init_states, u0 = self.shift_timestep(0.15, t0, init_states, u, self.f)

            state_history = np.dstack((state_history,self.DM2Arr(X0))) # 3D Matrix states values
            input_history = np.vstack((input_history,self.DM2Arr(u[:, 0]))) # 3D Input values
    
            t = np.vstack((t,t0))
            
            init_states = ca.reshape(init_states,-1,1)
            xx.append(init_states.full()) 

             ### LOS ####

            self.eta[0] = int(xx[0][0])
            self.eta[1] = int(xx[0][1])
            self.eta[5] = int(xx[0][2]) 

            X0 = ca.horzcat(X0[:, 1:],ca.reshape(X0[:, -1], -1, 1))


            t2 = time.time()

            times = np.vstack((
                times,
                t2-t1
            ))

            mpc_iteration=mpc_iteration+1
            print(mpc_iteration)

        ss_error = ca.norm_2(init_states - target_states)

        print(ss_error)

        # result = {'state_history':state_history,'init_history':input_history,'usv_state':xx,'t':t,'init_states':init_states,'target_states':target_states}
        
        output = dict(state_history=state_history,usv_states=xx,t=t,u_c=u_c,init_states=init_states,target_states=target_states,home_pose=home_pose,control_signal=u[:,0])
        return output

        # return state_history,xx,t,init_states,target_states,home_pose,u
    
        
if __name__ == "__main__":
    mpc = MPCUSV()
    # state_history,xx,t,init_states,target_states,home_pose,u=mpc.execute_MPC()
    output=mpc.execute_MPC()
    print('u_c',output['u_c'])

    print('u_c shape',np.shape(output['u_c']))

    
    Draw_MPC_point_stabilization(output['state_history'],output['usv_states'],output['home_pose'],output['target_states'])
    





