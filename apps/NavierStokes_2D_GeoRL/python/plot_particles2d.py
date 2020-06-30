import numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

import sys
if sys.version_info < (3, 5):
    raise "must use python 3.5 or greater"
import os
import os.path
from os import path

import subprocess

def execute(command):
    return subprocess.run([command, '-l'], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')


def readFromFile(name,newshape=None, typ = np.float32, numberOfItems = -1, seek_num = -1):
    ff = open(name, 'rb')
    if seek_num>0:
        ff.seek(seek_num)
    data = np.fromfile(ff, dtype=typ, count = numberOfItems) # numberOfItems = -1 means all items
    #print(name, data.shape, newshape)
    ff.close()
    if not(newshape == None):
        data.shape = newshape
    return data

def readConf(filename):
    prec=np.float64
    seek_num = 12
    llx, lly, urx, ury, nu, population_threshold, bc_treshold_xl, bc_treshold_xr, bc_treshold_yl, bc_treshold_yr, Uinfty_x, Uinfty_y = readFromFile(filename, None, prec, seek_num)

    if prec==np.float32:
        seek_num *= 4
    if prec==np.float64:
        seek_num *= 8

    remesh_ison, remesh_steps, remesh_method, N, nx, ny, order_ODEsolver, bc_xl, bc_xr, bc_yl, bc_yr = readFromFile(filename, None, np.int32, -1, seek_num)

    domain = np.array((llx, urx, lly, ury))
    remesh_ison = bool(remesh_ison)
    remesh_method = int(remesh_method)

    return domain, nx, ny, nu

def getFields(filename, timestep=0):

    filename+="_t"+str(timestep)

    domain, nx, ny, nu = readConf(filename+"_conf.dat")

    #if (N!= nx*ny):
    #    reshape = None
    #else:
    reshape = (ny, nx)

    prec=np.float64

    fn = filename+"_posx"+".dat"
    x = readFromFile(fn, reshape, prec)
    fn = filename+"_posy"+".dat"
    y = readFromFile(fn, reshape, prec)
    fn = filename+"_Ux.dat"
    Ux = readFromFile(fn, reshape, prec)
    fn = filename+"_Uy.dat"
    Uy = readFromFile(fn, reshape, prec)
    fn = filename+"_omega.dat"
    omega = readFromFile(fn, reshape, prec)
    fn = filename+"_time"+".dat"
    time = readFromFile(fn, None, prec)
    return domain, x, y, Ux, Uy, omega, nu, time, nx, ny

def getTime(filename, timestep):
    prec=np.float64
    fn=filename+"_t"+str(timestep)+"_time"+".dat"
    time = readFromFile(fn, None, prec)
    return time


def readParticles(filename, timestep, velocityField=False, backgroundField=False):
    prec=np.float64

    domain, nx, ny, N, nu, time = readConf(filename, timestep)

    if (N!= nx*ny):
        reshape = None
    else:
        reshape = (ny, nx)

    fn = filename+"_t"+str(timestep)+"_posx"+vel_str+".dat"
    x = readFromFile(fn, reshape, prec)
    fn = filename+"_t"+str(timestep)+"_posy"+vel_str+".dat"
    y = readFromFile(fn, reshape, prec)
    if velocityField:
        #reshape = None
        fn = filename+"_t"+str(timestep)+"_Ux.dat"
        Ux = readFromFile(fn, reshape, prec)
        fn = filename+"_t"+str(timestep)+"_Uy.dat"
        Uy = readFromFile(fn, reshape, prec)
        if backgroundField:
            fn = filename+"_t"+str(timestep)+"_Ux_bg.dat"
            Ux_bg = readFromFile(fn, reshape, prec)
            fn = filename+"_t"+str(timestep)+"_Uy_bg.dat"
            Uy_bg = readFromFile(fn, reshape, prec)
            return domain, x, y, Ux, Uy, Ux_bg, Uy_bg, nu, time, nx, ny
        else:
            return domain, x, y, Ux, Uy, nu, time, nx, ny
    else:
        fn = filename+"_t"+str(timestep)+"_omega.dat"
        omega = readFromFile(fn, reshape, prec)
        return domain, x, y, omega, nu, time, nx, ny

if __name__== "__main__":
    filename="test"
    command= "cd ../build/ && make "
    tmp=execute(command).strip()
    command= "../build/init --of "+str(filename)+" --nu 0.0 --domain -1 -1 1 1 --Uinfty 0.005 0 --nx 100 --ny 100"
    print("executing the command", command)
    tmp=execute(command).strip()

    iterations=20000
    final_time=100.
    for i in range(0,iterations):
        #print("================================== back in python, iteration number ",i)
        command= "../build/app --if test_t"+str(i)+" --of "+str(filename)+" --dt 0.1 --T "+str(final_time)+" --origo -.5 0 --semimajoraxis .25 --semiminoraxis .15"
        print("executing the command", command)
        tmp=execute(command).strip()
        time = getTime(filename,i+1)
        print("time=", time)
        if time>= final_time:
            iterations=i
            break

    domain, x, y, Ux, Uy, omega, nu, time, nx, ny = getFields(filename, iterations)
    im = plt.imshow(omega, extent=domain, origin="lower")#, vmin=-40, vmax=40)#, cmap=blue_red1)
    plt.colorbar()
    plt.savefig("test_omega.png")

    

