from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
import glob, math, os

ax = Axes3D(plt.figure())
os.chdir('3dp_simd_mpi')
norm = plt.Normalize()
nfiles = len(glob.glob("cells*.dat"))
colors = plt.cm.jet(norm(range(0,nfiles)))
colors[:,3] = 0.1
for file in glob.glob("cells*.dat"):
    data = np.genfromtxt(file, delimiter=' ')
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    r = data[:,3]
    numChilds = data[:,4]
    for i in range(0, len(x)):
        print file,i
        if numChilds[i] < 1e-6:
            X = np.array([[x[i]-r[i], y[i]-r[i], z[i]-r[i]],
                          [x[i]+r[i], y[i]-r[i], z[i]-r[i]],
                          [x[i]+r[i], y[i]+r[i], z[i]-r[i]],
                          [x[i]-r[i], y[i]+r[i], z[i]-r[i]],
                          [x[i]-r[i], y[i]-r[i], z[i]+r[i]],
                          [x[i]+r[i], y[i]-r[i], z[i]+r[i]],
                          [x[i]+r[i], y[i]+r[i], z[i]+r[i]],
                          [x[i]-r[i], y[i]+r[i], z[i]+r[i]]])
            face = [[X[0],X[1],X[2],X[3]],
                    [X[4],X[5],X[6],X[7]],
                    [X[0],X[1],X[5],X[4]],
                    [X[2],X[3],X[7],X[6]],
                    [X[1],X[2],X[6],X[5]],
                    [X[4],X[7],X[3],X[0]],
                    [X[2],X[3],X[7],X[6]]]
            ax.scatter3D(X[:,0], X[:,1], X[:,2], c=colors[int(file[5:9])])
            ax.add_collection3d(Poly3DCollection(face, facecolors=colors[int(file[5:9])]))
plt.show()
