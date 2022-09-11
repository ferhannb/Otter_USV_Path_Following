import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

df = pd.read_csv('otterRadiuesValues.csv')
df.drop(df[df['avarage_radius'] >= 1000].index, inplace = True)
# df.drop(df[df['iskele'] <= 0].index, inplace = True)
# df.drop(df[df['sancak'] <= 0].index, inplace = True)
iskele=df["iskele"]
sancak=df["sancak"]
radius=np.array([df["avarage_radius"].where(df['avarage_radius']<=100,100)])



X,Y = np.meshgrid(iskele,sancak)

print(X,Y)
ax = plt.axes(projection='3d')
# ax.axes.set_zlim3d(bottom=-500, top=500) 
ax.scatter(iskele, sancak, radius);
plt.show()


