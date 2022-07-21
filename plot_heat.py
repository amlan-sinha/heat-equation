from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import imageio
import os
import csv

plt.switch_backend('agg')

# Gifs
def generate_gif(filenames,output_path):
    with imageio.get_writer(output_path, mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

################################################################################################################################
############################################################# MAIN #############################################################
################################################################################################################################

path = './results/'
sysparams = np.loadtxt('./sysparams.txt',skiprows=1)

# Grid
xpts   = int(sysparams[0])
x_max  = float(sysparams[1])
x_min  = float(sysparams[2])
x_vec  = np.linspace(x_min,x_max,xpts+1)
dx     = (x_max-x_min)/xpts

ypts   = int(sysparams[3])
y_max  = float(sysparams[4])
y_min  = float(sysparams[5])
y_vec  = np.linspace(y_min,y_max,ypts+1)
dy     = (y_max-y_min)/ypts

numfiles = 0
for file in os.listdir(path):
        if file.endswith(".txt"):
            numfiles += 1

t_vec = np.zeros((numfiles))

for i in range(0,numfiles):
    with open(path+'t'+str(i)+'.txt') as fp:
        t = fp.readline()
        t_vec[i] = float(t)
        next(fp)
        csvfile = csv.reader(fp, dialect='excel-tab')
        for count,row in enumerate(csvfile):
            tot_len = count

tot_len += 1
x_vec = np.zeros((tot_len))
y_vec = np.zeros((tot_len))
U_vec = np.zeros((tot_len))
umax  = 0

datapath = path+'U(t,x,y)'

try:
    os.makedirs(datapath)
except OSError:
    print('Was unable to create new directory, it may already exist.')
else:
    print('Was able to create new directory')

for i in range(0,numfiles):
    
    data = np.loadtxt(path+'t'+str(i)+'.txt',skiprows=2)
    x_vec = data[0:-1,0]
    y_vec = data[0:-1:,1]
    U_vec = data[0:-1:,2]

    umax = np.amax(U_vec) if np.amax(U_vec)>umax else umax;

    print('Current maximum: %f' %umax)

for i in range(0,numfiles):

    print('Currently processing file %d' %i)

    data = np.loadtxt(path+'t'+str(i)+'.txt',skiprows=2)
    x_vec = data[0:-1,0]
    y_vec = data[0:-1:,1]
    U_vec = data[0:-1:,2]

    fig  = plt.figure(1)
    ax   = fig.add_subplot(1,1,1, projection='3d')
    surf = ax.plot_trisurf(x_vec, y_vec, U_vec, cmap=cm.viridis, antialiased=False, vmin=0.0, vmax=umax)
    fig.colorbar(surf)
    ax.set_title('time = %3.3f s' % t_vec[i])
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$U(t,x,y)$')
    ax.set_zlim(0.0,umax)
    plt.savefig(path+'U(t,x,y)/t'+str(i).zfill(3)+'.png', dpi=300)
    plt.clf()

filenames = [path+'U(t,x,y)/t'+str(i).zfill(3)+'.png' for i in range(0,numfiles)]
generate_gif(filenames,path+'U(t,x,y)/U(t,x,y).gif')
