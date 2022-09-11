
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import pchip_interpolate



class spline_gen():

        def __init__(self):
                self.Wpx = []
                self.Wpy = []

        def pchip_gen(self,Wpx,Wpy):
                self.Wpx=Wpx
                self.Wpy=Wpy

                min_Wpx = min(Wpx,key=lambda x: x) # min waypoint değeri
                max_Wpx = max(Wpx,key=lambda x: x) # max waypoint değeri

                self.x_interpolate = np.linspace(min_Wpx,max_Wpx,100*(len(Wpx)-1))
                self.y_interpolate = pchip_interpolate([wpx for wpx in Wpx],[wpy for wpy in Wpy],self.x_interpolate)
                print(np.shape(self.y_interpolate))

                return self.x_interpolate,self.y_interpolate

        def f(self,x,f):
                return f[0]*x**3+f[1]*x**2+f[2]*x+f[3] 


        def execute(self):
                m=0
                self.c_list=[]
                index_list=[]
                for i in range(len(self.Wpx)-1):
                        indx=int((self.Wpx[i+1]-self.Wpx[i])/(abs(self.Wpx[0]-self.Wpx[-1]))*len(self.y_interpolate))
                        index_list.append(indx)
                        print(indx)
                        x = np.linspace(self.Wpx[i],self.Wpx[i+1],index_list[i])
                        if i ==0:

                                y = self.y_interpolate[0:indx]
                        else:
                                m = np.sum(index_list)
                                n = np.sum(index_list[0:(i)])
                                y = self.y_interpolate[n:m]
                        print(index_list)
                        coeff = np.polyfit(x,y,3)
                        self.c_list.append(coeff)
                        plt.plot(self.Wpx,self.Wpy,'r.')
                        plt.plot(x,self.f(x,coeff),'g--')
                        
                        plt.show()

                print(type(self.c_list))
                return self.c_list 

                
if __name__=='__main__':


        Wpx=[4,12,12.5,34,40]
        Wpy=[4,20,25,67,45]


        los=spline_gen()
        xint,yint =los.pchip_gen(Wpx,Wpy)
        plt.plot(Wpx,Wpy,'ro')
        plt.plot(xint,yint,'b.')
        plt.show()
        los.execute()
