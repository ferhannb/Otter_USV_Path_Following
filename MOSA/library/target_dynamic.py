
import numpy as np
import math
import matplotlib.pyplot as plt

class Target():

    def __init__(self,X=0,Y=0,theta=45,V=2,omega=0,dt=0.1):

        self.V=V
        self.omega=omega
        self.x=X
        self.y=Y
        self.theta=theta
        self.dt=dt

    def movement(self):

        A = np.array([[math.cos(math.radians(self.theta)),0],
                      [math.sin(math.radians(self.theta)),0],
                      [0,1]])

        B = [[self.V],[self.omega]]

        X_ = [[self.x],[self.y],[self.theta]]

        X_ = X_ + self.dt*(A @ B)

        return X_



if __name__=='__main__':

    target = Target()
    X = [[0],[0],[0]]
    x = []
    y = []
    for i in range(200):
        print(i)
        X = X+target.movement()
        x.append(X[0])
        y.append(X[1])

    print(X)
    plt.plot(x,y)
    plt.show()


    




