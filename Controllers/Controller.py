

import casadi as ca 
import numpy as np
from OtterDynamicClass import Otter



class MPCController(Otter):

    def __init__(self, N=10, dt = 0.1,Q=[1,0,0,10,10],R=[1,1]):
        Otter.__init__(self)
        self.predict_st = np.array([[0],[0],[0],[0],[0]])
        self.N = N      # Horizon
        self.dt = dt    # Sampling Time 
         
        self.prediction_st=np.array([0,0,0,0,0])
        self.real_states = np.array([[self.nu[0]], [self.nu[5]],[self.current_eta[0]],[self.current_eta[1]],[self.current_eta[5]]])    
        self.u      = ca.SX.sym('u')
        self.r      = ca.SX.sym('r')
        self.x      = ca.SX.sym('x')
        self.y      = ca.SX.sym('y')
        self.psi    = ca.SX.sym('psi')
        self.states =ca.vertcat(self.u,self.r,self.x,self.y,self.psi)  # state vector
        self.n_states =self.states.numel()  # number of states

        # inputs otter
        self.Tp = ca.SX.sym('Tp')  # Port input (Newtoon)
        self.Ts = ca.SX.sym('Ts')  # Starboard input (Newtoon)
        
        self.inputs = ca.vertcat(self.Tp,self.Ts) # input vector
        self.B = 1.08/2 # Beam  
        self.n_inputs = self.inputs.numel() # number of inputs

        ### Dynamic equation of Otter ###
        self.dynamics_otter = ca.vertcat( (-45.3825*self.u+1.7315*ca.fabs(self.u)*self.u+0.5257*(self.Tp+self.Ts))/(50-12.776),
                                           -1.48517*self.r-0.8443*self.r*ca.fabs(self.r)-10.7*self.r*ca.fabs(self.r)+0.019*self.B*(self.Tp-self.Ts),
                                           self.u*ca.cos(self.psi),
                                           self.u*ca.sin(self.psi),
                                           self.r)
        
        ## STATE WEIGHT MATRIX ##
        self.Q = ca.SX.zeros(self.n_states,self.n_states) 
        self.Q[0,0] = Q[0]
        self.Q[1,1] = Q[1]
        self.Q[2,2] = Q[2]
        self.Q[3,3] = Q[3]
        self.Q[4,4] = Q[4]
        
        self.Q_last = ca.SX.zeros(self.n_states,self.n_states) 
        self.Q_last[0,0] = Q[0]#2.2
        self.Q_last[1,1] = 0#2.2
        self.Q_last[2,2] = Q[2]#1
        self.Q_last[3,3] = Q[3]#1
        self.Q_last[4,4] = Q[4]#1

        ## INPUT WEIGHT MATRIX
        self.R = ca.SX.zeros(self.n_inputs,self.n_inputs) 
        self.R[0,0] =R[0]
        self.R[1,1] =R[1]

        self.f = ca.Function('f',[self.states,self.inputs],[self.dynamics_otter]) 
        

        self.X = ca.SX.sym('X',self.n_states,self.N+1)        # Vector representing all states in time up to the horizon to be controlled.
        self.U = ca.SX.sym('U',self.n_inputs,self.N)          # Vector representing all inputs in time up to the horizon to be controlled.
        self.P = ca.SX.sym('P',self.n_states+self.n_states)   # Vector representing target and initial states.

        self.cost_function = 0
        self.g = self.X[:, 0] - self.P[:self.n_states] 
        ## Runge-Kutta opt.
        for k in range(self.N):
            st = self.X[:, k]
            input = self.U[:, k]
            if k < self.N-2:
                delta_u = self.U[:,k+1]-input
            if k<=self.N-2:
                self.cost_function = self.cost_function + (st - self.P[self.n_states:]).T @ self.Q @ (st - self.P[self.n_states:]) + delta_u.T @ self.R @ delta_u
            if k==self.N-1:
                self.cost_function = self.cost_function + (st- self.P[self.n_states:]).T @ self.Q_last @ (st- self.P[self.n_states:])   # Terminal Term 
            st_next = self.X[:, k+1]
            k1 = self.f(st, input)
            k2 = self.f(st + self.dt/2*k1, input)
            k3 = self.f(st + self.dt/2*k2, input)
            k4 = self.f(st + self.dt * k3, input)
            st_next_RK4 = st + (self.dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)- (self.real_states - self.predict_st)
            st_next_RK4 = st_next_RK4 
            self.g = ca.vertcat(self.g, st_next - st_next_RK4)
        


            #######################################

            ## Constraints 
            
        

        OPT_variables = ca.vertcat( ca.reshape(self.X, -1, 1), ca.reshape(self.U, -1, 1))

        self.nlp_prob = {
                    'f': self.cost_function,
                    'x': OPT_variables,
                    'g': self.g,
                    'p': self.P
                    }


        opts = {
                'ipopt': {
                    'max_iter': 500,
                    'print_level': 0,
                    'acceptable_tol': 1e-8,
                    'acceptable_obj_change_tol': 1e-6
                },
                'print_time': 0
            }
            
        ## State ve İnput değerlerinin Kısıtlandırılması 
        self.solver = ca.nlpsol('solver','ipopt',self.nlp_prob,opts)


        lbx = ca.DM.zeros((self.n_states*(self.N+1) + self.n_inputs*self.N, 1))
        ubx = ca.DM.zeros((self.n_states*(self.N+1) + self.n_inputs*self.N, 1))

        lbx[0: self.n_states*(N+1): self.n_states] = -1.0          # u lower bound
        lbx[1: self.n_states*(N+1): self.n_states] = -0.35         # r lower bound
        lbx[2: self.n_states*(N+1): self.n_states] = -ca.inf       # x lower bound
        lbx[3: self.n_states*(N+1): self.n_states] = -ca.inf       # y lower bound
        lbx[4: self.n_states*(N+1): self.n_states] = -ca.inf       # psi lower bound

        ubx[0: self.n_states*(N+1): self.n_states] = 1            # u lower bound
        ubx[1: self.n_states*(N+1): self.n_states] = 0.35         # r lower bound
        ubx[2: self.n_states*(N+1): self.n_states] = ca.inf       # x lower bound
        ubx[3: self.n_states*(N+1): self.n_states] = ca.inf       # y lower bound
        ubx[4: self.n_states*(N+1): self.n_states] = ca.inf       # psi lower bound

 
        lbx[self.n_states*(self.N+1)::self.n_inputs] =   -70     # Tp lower bound
        ubx[self.n_states*(self.N+1)::self.n_inputs] =    120    # Tp upper bound
        lbx[self.n_states*(self.N+1)+1::self.n_inputs] = -70     # Ts lower bound
        ubx[self.n_states*(self.N+1)+1::self.n_inputs] =  120    # Tp upper bound
        
        
        
        # equality constraints 

        lbg = ca.DM.zeros(self.n_states*(self.N+1),1)
        ubg = ca.DM.zeros(self.n_states*(self.N+1),1)
       

        self.args = {
                'lbg': lbg,  # constraints lower bound
                'ubg': ubg,  # constraints upper bound
                'lbx': lbx,
                'ubx': ubx
                    }
    
        
    def execute_MPC(self,cnt = 0, init_states=[0,0,0,0,0],target_states=[3,0,0,0,0]):

        init_states = np.array(init_states)
        self.heading_error = abs(init_states[2]-target_states[2])
        u0 = ca.DM.zeros((self.N,self.n_inputs))
        X0 = ca.repmat(init_states,1,self.N+1)
                     
        u_control_prediction=[] 

        self.args['p'] = ca.vertcat(
            init_states,    # current state
            target_states   # target state
        )
    
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
        u_control_prediction.append(u[:,0:self.N])
        u_control_prediction = np.array(u_control_prediction)
        u_control_prediction.shape = (2,self.N)

        X0 = ca.reshape(sol['x'][: self.n_states * (self.N+1)], self.n_states, self.N+1)
        u_contol = np.array(u[:,0].full())
        u_contol.shape=(2,)
        self.prediction_st = np.array(X0[:,0].full())
        self.predict_st = np.array(X0[:,0].full())
        X0 = ca.horzcat(X0[:, 1:],ca.reshape(X0[:, -1], -1, 1))
        self.u_control = u_contol
        self.function()

      
        output = dict(state_prediction=self.prediction_st[0],input_prediction = u_control_prediction, control_signal=u_contol,target_states=target_states)
        return output