from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import glob, math, os

ax = Axes3D(plt.figure())
os.chdir('tests')
norm = plt.Normalize()
nfiles = len(glob.glob("bodies*.dat"))
colors = plt.cm.jet(norm(range(0,nfiles)))
for file in glob.glob("bodies*.dat"):
    data = np.genfromtxt(file, delimiter=' ')
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    ran = range(0, len(x), 1)
    ax.scatter(x[ran], y[ran], z[ran], c=colors[int(file[6:10])])
plt.show()
