import numpy as np 
import matplotlib.pyplot as plt


class FossenCubicSpline():

    def __init__(self,x,y):

        self.Wpx = x
        self.Wpy = y


    def Spline_generator(self,xstart,xfinal,ystart,yfinal):

        N = len(self.Wpx)-1 # Toplam Denklem Sayısı
        print(N)
        th = list(range(0,(N+1)))  # theta values along the path
       
        c1 = [3*th[0]**2,2*th[0],1,0]

        c2 = [3*th[N]**2,2*th[N],1,0]
        print('c2',c2)
        ## Y vector generation accorting  to number of waypoint

        # xstart = self.Wpx[0]
        # xfinal = self.Wpx[-1]

        # ystart = self.Wpy[0]
        # yfinal = self.Wpy[-1]


        # y_x = np.zeros((4*N,1))
        # y_y = np.zeros((4*N,1))

        y_x = [xstart,self.Wpx[0]]
        y_y = [ystart,self.Wpy[0]]

        ## Generate A Matrix

        for k in range(1,N):
            
            y_x.extend([self.Wpx[k],self.Wpx[k],0,0])
            y_y.extend([self.Wpy[k],self.Wpy[k],0,0])
        
        y_x.extend([self.Wpx[N],xfinal])
        y_y.extend([self.Wpy[N],yfinal])

        A = np.zeros((4*N,4*N))
        A[0,0:4] = c1
        A[1,0:4] = np.array([th[0]**3,th[0]**2,th[0],1])
        A[4*N-1,(4*(N-1)):4*N] = c2
        A[(4*N-1)-1,(4*(N-1)):4*N] = [th[N]**3,th[N]**2,th[N],1]


        for j in range(1,N):

            row = 4*(j-1)+2
            col = 4*(j)
            print(np.shape( A[row+2,col-4:col]))
            A[row,col-4:col] = [th[j]**3,th[j]**2,th[j],1]
            A[row+1,col:col+4] = [th[j]**3,th[j]**2,th[j],1]
            
            A[row+2,col-4:col] = [-3*th[j]**2,-2*th[j],-1,0]
            A[row+2,col:col+4] = [3*th[j]**2,2*th[j],1,0]
            
            A[row+3,col-4:col] = [-6*th[j],-2,0,0]
            A[row+3,col:col+4] = [6*th[j],2,0,0]

        x_x = np.matmul(np.linalg.inv(A),y_x)
        
        x_y = np.matmul(np.linalg.inv(A),y_y)

        Acoeff = np.zeros((N,N))
        Bcoeff = np.zeros((N,N))

        for m in range(0,N):

            Acoeff[m,:] = x_x[4*m:4*m+4]
            Bcoeff[m,:] = x_y[4*m:4*m+4]

            # x_theta_one[m,:] = [0,3*Acoeff[m,0],2*Acoeff[m,1],Acoeff[m,2]]
            # y_theta_one[m,:] = [0,3*Bcoeff[m,0],2*Bcoeff[m,1],Bcoeff[m,2]]

            # x_theta_two[m,:] = [0,0,6*Acoeff[m,0],2*Acoeff[m,1]]
            # y_theta_two[m,:] = [0,0,6*Bcoeff[m,0],2*Bcoeff[m,1]]

            # x_theta_three[m,:] = [0,0,0,6*Acoeff[m,0]]
            # y_theta_three[m,:] = [0,0,0,Bcoeff[m,0]]


    ## Path Length Calculation
        path_length = 0
        for i in range(0,N):

            th = np.linspace(i,i+1,100)
            
            xpath = np.polyval(Acoeff[i,:],th)
            ypath = np.polyval(Bcoeff[i,:],th)

            for k in range(len(xpath)-1):
                path_length += np.sqrt((xpath[k+1]-xpath[k])**2+(ypath[k+1]-ypath[k])**2) 

     
        print(path_length)



        




        
if __name__ == '__main__':
    x = [0, 200, 400, 700 ,1000]
    y = [0, 600, 500, 400, 1200]
    cubic = FossenCubicSpline(x,y)
    cubic.Spline_generator(0,0,0,0)




