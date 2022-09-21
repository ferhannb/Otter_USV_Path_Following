#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import clip
import math

class USVController():
    def __init__(self):
        self.prev_error = 0
        self.error = 0
        self.prev_heading_error = 0
        self.prev_speed_error = 0
        self.delta_error = 0
        self.saturation_limit = 0
        self.U_i_heading = 0
        self.U_i_speed = 0
        self.U_i_old = 0
        self.filtred_signal = 0.1
        

    def Heading_controller(self,ref_heading,current_eta,Saturation_limit=104,Kp=0,Ki=0,Kd=0):

        
        feed_back = math.degrees(current_eta[-1])
        heading_error = ref_heading-feed_back

        while heading_error<-180:
            heading_error = heading_error +360
        while heading_error > 180:
            heading_error = heading_error -360

        # Proportion Term
        U_p = Kp*heading_error

        # Integral Term
        if self.U_i_heading < Saturation_limit:
            self.U_i_heading += Ki*heading_error 
        else:
            pass

        # Derivative Term
        U_d = Kd*(heading_error-self.prev_heading_error)

        

        # Total Control Signal 
        Control_signal = U_p + self.U_i_heading + U_d 

        self.prev_heading_error  = heading_error

        return Control_signal 


        
    def Speed_controller(self,velocity,ref_speed,Saturation_limit=104,Kp=0,Ki=0,Kd=0,Kf=0):
        Umin = -104
        Umax = 104
        speed_error = velocity-ref_speed
        # Feedforward Term
        U_f = 0

        # Proportion Term 
        U_p = Kp*speed_error
        # Derivative Term 
        U_d = Kd*(speed_error-self.prev_speed_error)
        self.prev_speed_error = speed_error
        # Integral Term 
        self.U_i_speed =self.U_i_old + Ki*speed_error
        self.U_i_old = self.U_i_speed
        # Control Signal 
        U = U_p + U_d + self.U_i_speed
        # Wind-up
        U_filtred = clip(U,Umin,Umax)
        if U!=U_filtred:
             self.U_i_old = U_filtred-U_p-U_f

        return U_filtred


    def control_allocation(self,u_avr,u_diff):
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


    def Referans_speed_signal_filter(self,speed_ref):

        delta_rate = 0.03
        if self.pre_speed_ref == self.filtred_signal:
            self.pre_speed_ref = self.filtred_signal
        else:
            if self.filtred_signal<speed_ref:
                self.filtred_signal += +delta_rate
                if self.filtred_signal>speed_ref:
                    self.filtred_signal=speed_ref

            elif self.filtred_signal>speed_ref:
                self.filtred_signal += - delta_rate 
                if self.filtred_signal < speed_ref:
                    self.filtred_signal = speed_ref
            else:
                self.filtred_signal=speed_ref
        return self.filtred_signal


    def Filtred_heading_referans(heading_ref,filtred_heading_signal):

        if heading_ref-filtred_heading_signal<0:
            if (heading_ref-filtred_heading_signal)%360<(filtred_heading_signal-heading_ref):
                filtred_heading_signal = filtred_heading_signal+0.36

                if filtred_heading_signal>0:

                    if (-filtred_heading_signal)%360>heading_ref:
                        filtred_heading_signal=heading_ref
            elif  (heading_ref-filtred_heading_signal)%360>(filtred_heading_signal-heading_ref):

                filtred_heading_signal = filtred_heading_signal - 0.36
                if filtred_heading_signal<0:
                    if (filtred_heading_signal%360)>heading_ref:
                        filtred_heading_signal=heading_ref
        elif heading_ref-filtred_heading_signal>0:
            if (heading_ref-filtred_heading_signal)<(filtred_heading_signal-heading_ref)%360:
                filtred_heading_signal = filtred_heading_signal+0.36
                if (filtred_heading_signal)>0:
                    if filtred_heading_signal>heading_ref:
                        filtred_heading_signal=heading_ref
            elif  (heading_ref-filtred_heading_signal)>(filtred_heading_signal-heading_ref)%360:

                filtred_heading_signal = filtred_heading_signal - 0.36
                if filtred_heading_signal<0:

                    if (filtred_heading_signal%360)<heading_ref:
                            filtred_heading_signal=heading_ref
                
        return filtred_heading_signal


    
            



         
         


        

