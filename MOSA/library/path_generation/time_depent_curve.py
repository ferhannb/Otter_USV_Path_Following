

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import pchip_interpolate,Rbf,splrep,splev
from scipy.optimize import curve_fit


class spline_gen():

        def __init__(self):
                self.Wpx = []
                self.Wpy = []
                self.time = []

        def pchip_gen(self,Wpx,Wpy,Time):

                self.Wpx=Wpx
                self.Wpy=Wpy
                self.time=Time
                min_time = min(Time,key=lambda x: x) # min waypoint değeri
                max_time = max(Time,key=lambda x: x) # max waypoint değeri

                # X konum bilgileri
                self.time_interpolate = np.linspace(min_time,max_time,1000*(len(self.time)-1))
                self.x_interpolate = pchip_interpolate([time for time in Time],[wpx for wpx in Wpx],self.time_interpolate)
                
                # Y konum bilgileri
                self.time_interpolate = np.linspace(min_time,max_time,1000*(len(self.time)-1))
                self.y_interpolate = pchip_interpolate([time for time in Time],[wpy for wpy in Wpy],self.time_interpolate)

                return self.x_interpolate,self.y_interpolate
        


        def df(self,f):
            output = []
            for i in range(len(self.time_interpolate)):

                output=3*f[0]*self.time_interpolate[i]**2+2*f[1]*self.time_interpolate[i]+f[2]
                output.append(output)
            return output

        def ddf(self,f):
            return 6*f[0]*self.time_interpolate+2*f[1]

        def execute(self):
                m=0
                self.coef_x_list=[]
                self.coef_y_list=[]
                index_list_x=[]
                index_list_y=[]
                x_value=[]
                y_value=[]

                for i in range(len(self.Wpx)-1):
                        indx_x=int((self.time[i+1]-self.time[i])/(abs(self.time[0]-self.time[-1]))*len(self.x_interpolate))
                        indx_y=int((self.time[i+1]-self.time[i])/(abs(self.time[0]-self.time[-1]))*len(self.y_interpolate))
                        
                        indx_time = int(int((self.time[i+1]-self.time[i])/(abs(self.time[0]-self.time[-1]))*len(self.time_interpolate)))
                        index_list_x.append(indx_x)
                        index_list_y.append(indx_y)
                        x = np.linspace(self.time[i],self.time[i+1],index_list_x[i])
                        y = np.linspace(self.time[i],self.time[i+1],index_list_y[i])

                        if i ==0:

                                y_x = self.y_interpolate[0:indx_x]
                                y_y = self.y_interpolate[0:indx_y]
                        else:
                                m_x = np.sum(index_list_x)
                                n_x = np.sum(index_list_x[0:(i)])
                                y_x = self.y_interpolate[n_x:m_x]

                                m_y = np.sum(index_list_y)
                                n_y = np.sum(index_list_y[0:(i)])
                                y_y = self.y_interpolate[n_y:m_y]
                                print('****')

                        
                        coeff_x = np.polyfit(x,y_x,3)
                        coeff_y = np.polyfit(y,y_y,3)
                        self.coef_x_list.append(coeff_x)
                        self.coef_x_list.append(coeff_y)
                        
        
                return self.coef_x_list,self.coef_x_list



                
if __name__=='__main__':


        Wpx=[4,12,20,15,40,65,10]
        Wpy=[4,20,25,67,45,78,95]
        Time =[0,10,57,58,450,1110,1200]

        los=spline_gen()
        xint,yint =los.pchip_gen(Wpx,Wpy,Time)
        coef_x, coef_y = los.execute()

        coef_x_first=coef_x[0]

        plt.plot(Wpx,Wpy,'ro')
        plt.plot(xint,yint,'b-.')
        plt.show()
        los.execute()
