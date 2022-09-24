#!/usr/bin/env python
'''
    :project:   PID Control Library for Dynamic Systems
    :author:    email:manyetar@gmail.com
                github:saltin0

    :reference: Fossen T, Handbook of Marinecraft Systems and Motion Control
                Ogata K, Classical Control Theory
                Nisse N, Modern Control
'''
from math import sqrt

class PID():
    def __init__(self,Kp,Ki,Kd,Kp1=None,
                derivative_feedback=False,derivative_filter_coeff=None,
                inner_loop_p = False,
                saturation_limit=None,
                wind_feed_forward = False,
                Kf = None,
                acceleration_gain = None,
                wave_current_balancer = False):
        '''
            Initiate the PID controller structure parameters

            :params:    Kp-Proportional gain
            :params:    Ki-Integral gain
            :params:    Kd-Derivative gain
            :params:    derivative_feedback-Set the derivative gain in measurement way
                        PI-D Controller
            :params:    saturation_limit-Set the maximum value of controller corresponding actuator
            :params:    derivative_filter_coefficient-To set the system as causal system
            :params:    inner_loop_p-Change the system controller into the PIP controller. It is generally used
                        for first order systems.
            :params:    wind_feed_forward-Add the estimated wind disturbance to the controller. Estimation
                        must be done outside of the controller.
            :params:    Kf-Set the feedforward gain. Multiply the reference value and add it to the controller.
            :params:    acceleration_gain-To reduce the disturbance effect set the acceleration gain.
            :params:    wave_current_balance-To reduce the wave effect on the system set the variable as True with 
                        identified wave and RAO models.
            --------------------------------------------------------------------------------------------------------
            :variable:  system_output_pre-Keep the previous response of the system to take derivative
            :variable:  ref_pre- Keep the previous reference value (for derivative)
            :variable:  deriv_out_pre-Keep the previous derivative output of PID controller to reduce the derivative
                        gain with filter.
            :variable:  ref_now-Reference value holder.
            :variable:  err_cum-Cumulative error of the integral part of the PID.
            :variable:  clamp-Reset the integral at the boundaries via clamping method. Stop integrating 
                        at boundaries.
            :variable:  err_now-Ref_now - system_output_now
            :variable:  u-Keep the control signal in memory.
        '''
        # Parameters
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.Kp1 = Kp1
        self.derivative_feedback = derivative_feedback
        self.derivative_filter_coefficient = derivative_filter_coeff
        self.saturation_limit = saturation_limit
        self.Ag  = acceleration_gain 
        self.wind_ff = wind_feed_forward
        self.Kf = Kf
        self.inner_loop_p = inner_loop_p
        self.wave_current_balancer = wave_current_balancer
        # Variables
        self.system_output_pre = 0.0
        self.deriv_out_pre = 0.0
        self.ref_pre = 0.0
        self.ref_now = 0.0
        self.err_cum = 0.0
        self.clamp = False
        self.err_now = 0.0
        self.u = 0.0
        self.deriv_filter_integral_cum = 0.0

    def execute(self,system_output,dt,acceleration=0.0,nonlinear_terms=None):
        '''
            :function:  execute-Calculate the output of formed structure which obtains
                        feedforward-PID-extra P-nonlinear elimination-wind feed forward-
                        acceleration feedback - current compansator.
            
            :params:    system_output-Output of the plant model. Measured data.
            :params:    dt-Sample time.
            :params:    acceleration-Acceleration measurement. Data type is float.
            :params:    nonlinear_terms-Add nonlinearities of the system to eliminate.

            :return:    u-Output of the entire controller.
        '''
        self.system_output = system_output
        # Calculate error definitions 
        err_now = self.ref_now - system_output
        err_pre = self.ref_pre - self.system_output_pre

        ### Calculate parallel PID output ###
        proportion_part = self.propotional(err_now)
        integral_part = self.integrate(err_pre,err_now,dt)

        # Form the PI-D or PID blocks
        if self.derivative_feedback == True:
            derivative_part = self.derivate(self.system_output_pre,system_output,dt)
        if self.derivative_feedback == False:
            derivative_part = -1*float(self.derivate(err_pre,err_now,dt))
        # Sum the parallel structure
        u = float(integral_part + proportion_part + derivative_part)

        # If the system is controlled by PI-P1 structure then following code will be running
        if self.inner_loop_p:
            # Inner loop difference defines the difference between system output and
            #       PID controller output.
            inner_loop_diff = u - system_output
            u = self.in_proportional(inner_loop_diff)
        
        # Acceleration feedback mechanism added 
        #   this property will reduces the disturbance effect.
        if self.Ag and acceleration is not None:
            try:
                u = u - self.Ag*float(acceleration)
            except:
                pass
        # Add feedforward to control action
        if self.Kf is not None:
            u_kf =self.Kf*self.ref_now
            u = u+u_kf
        # Eliminate the nonlinearities of the system 
        if nonlinear_terms is not None:
            try:
                u = u - float(nonlinear_terms)
            except:
                pass

        # Anti-windup for integral part 
        #   Anti-windup method : clamping
        u_sat = self.reset_integral(u)
        
        if u_sat!=u:

            integral_part=u_sat-proportion_part-u_kf
        	
        ### End of PID ###

        # Return the parameters to access from outside for debug
        self.err_now = err_now
        
        # Recursively state keeping
        self.ref_pre = self.ref_now
        self.system_output_pre = system_output

        return u_sat

    def set_setpoint(self,reference):
        '''
            :function:  This function sets the reference signal and send to the controller
            :params:    reference-Set the reference signal controller may be achieve
        '''
        self.ref_now = reference

    def integrate(self,err_pre,err,dt):
        '''
            :function:  integrate
            :params:    err_pre-Error value one step before
            :params:    err-Error value for now

            :return:    err_cum-Cumulative error with multiplied by Ki

            Note :      Integration rule - Trapezoidal Rule
        '''
        # Check if the control signal saturated
        if self.clamp == False:
            # If saturation does not exist integrate the error term
            self.err_cum += self.Ki*(err_pre+err)*dt/2
        elif self.clamp == True:
            # If saturation exists do not integrate the error term 
            pass
        return self.err_cum
    
    def propotional(self,err):
        '''
            :function:  proportional
                        This function multiplies error by Kp
            
            :params:    err-Error for now

            :return:    self.Kp*err-Gained error value
        '''
        return self.Kp*err
    
    def in_proportional(self,err):
        '''
            :function:  proportional
                        This function multiplies error by Kp
            
            :params:    err-Difference between control signal and
                        system output.

            :return:    self.Kp*err-Gained error value
        '''
        return self.Kp1*err


    def derivate(self,pre_val,val,dt):
        '''
            :function:  This function takes the derivative of error or system output

            :params:    pre_val-Previous state for numerical derivation.
            :params:    val-Current value for derivation
            :params:    dt-Sampling time.

            :return:    Filtered derivative value.
        '''
        if self.derivative_filter_coefficient is not None:
            deriv_out = (self.Kd*(val-pre_val)-self.deriv_filter_integral_cum)*self.derivative_filter_coefficient
            self.deriv_filter_integral_cum += dt*deriv_out
            print("Deriv filter int : {}".format(self.deriv_filter_integral_cum))

            if dt > 0.0:
                return deriv_out
            else:
                return 0.0
        else:
            if dt > 0.0:
                return self.Kd*(val-pre_val)/dt
            else:
                return 0.0

    def reset_integral(self,u):
        '''
            :function:  reset_integral
                        This function saturates the controller output.
                        When the controller output is higher or lower than limits
                        clamping flag rises and stops the integrating process.
            
            :params:    u-Controller signal
            :return:    u-Saturated signal
        '''
        # Saturate the limit
        if self.saturation_limit is not None and u>self.saturation_limit:
            u = self.saturation_limit
            # If saturation is observed, rise clamp flag
            self.clamp = True
        elif self.saturation_limit is not None and u<-self.saturation_limit:
            u = -self.saturation_limit
            self.clamp = True
        else:
            self.clamp = False
        self.u = u
        return u

    '''
        Ocean current balancer.
    '''
    def wave_current_to_amplitude(self):
        pass

    def response_amplitude_operator(self):
        '''
            :function:      response_amplitude_operator: RAO operator to estimate 
                            the wave forces impact to the USV according to the wave
                            amplitude. Controller has to compansate the 
            
            :return:        wave_force_est: Estimated wave forces.
        '''
        wave_force_est = None
        return wave_force_est
    @property
    def see_err(self):
        return self.err_now

    @property
    def debugger(self):
        '''
            :property:  debugger-See all current and memorized data in controller.
        '''
        debug_vals = {
            "clamp"                 :   self.clamp,
            "err_now"               :   self.err_now,
            "err_cum"               :   self.err_cum,
            "ctrl_output"           :   self.u,
            "acc_feedback_gain"     :   self.Ag,
            #"saturation_limit"      :   self.saturation_limit,
            #"reference_now"         :   self.ref_now,
            #"reference_pre"         :   self.ref_pre,
            #"sys_out_pre"           :   self.system_output_pre,
            #"sys_out_now"           :   self.system_output,
            "Kp,Ki,Kd,Kp1"          :   (self.Kp,self.Ki,self.Kd,self.Kp1),
            "derivative_filter_coef":   self.derivative_filter_coefficient,
        }
        return debug_vals
    @property
    def get_lib_usage(self):
        lib_usage = "Explanation of PID library for USV's \
            First you need to construct the class in your code. Constructor takes the parameters as: \
            Kp,Ki,Kd as basic PID parameters and Kp1 for PIDP structure. In PIDP you first close \
            the loop with just P type controller to make the systm faster and close the second loop with PID \
            to force the system to act as desired. In derivative feedbak structure D type controller is \
            on the negative feedback way. To reduce the derivative gain set derivative feedback coefficient \
            this coefficient reduces the gain in high frequencies. Inner-loop P controller exist ot not parameters will be set \
            default is false. Set the saturation limit the save the actuators life and give the inputs on limits \
            This parameters automaticly set the anti-windup structure as clamping. If there is a wind sensor \
            you can calculate the wind effect in units of moment or forces and add into the system. Kf is feedforward \
            gain. Set the gain as de"
        return lib_usage

