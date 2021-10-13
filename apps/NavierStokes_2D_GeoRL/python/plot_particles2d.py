import numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap

import sys
if sys.version_info < (3, 5):
    raise "must use python 3.5 or greater"
import os
import os.path
from os import path

import subprocess

cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 1.0),
                   (0.5, 0.1, 0.0),
                   (1.0, 0.0, 0.0))
        }
blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)


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

if __name__== "__init__":
    build_dir = "../../../build"
    app_name = os.path.join(build_dir, "apps/NavierStokes_2D_GeoRL/navier_stokes_2d_georl")
    init_name = os.path.join(build_dir, "apps/NavierStokes_2D_GeoRL/init")

    ###
    # Reynolds number is defined as
    # Re = u L / nu
    # For a Reynolds number of 100 this means nu = 0.1*.15/100
    Re=500
    ux=0.1
    radius=0.15
    nu=ux*(2*radius)/Re
    ###
    filename="Re500_100x50"
    command= f"cd {build_dir} && make "
    tmp=execute(command).strip()
    command= f"{init_name} --of "+filename+" --nu "+str(nu)+" --domain -1 -1 3 1 --Uinfty "+str(ux)+" 0 --nx 100 --ny 50 --order 1"
    # --population_threshold 0.001"
    print("executing the command", command)
    tmp=execute(command).strip()

    iterations=20000
    final_time=50.
    for i in range(0,iterations):
        #print("================================== back in python, iteration number ",i)
        command= f"{app_name} --if "+filename+"_t"+str(i)+" --of "+filename+" --T "+str(final_time)+" --origo -.5 0.025 --semimajoraxis "+str(radius)+" --semiminoraxis "+str(radius)
        print("executing the command", command)
        tmp=execute(command).strip()
        time = getTime(filename,i+1)
        print("time=", time)
        if time>= final_time:
            iterations=i
            break

        domain, x, y, Ux, Uy, omega, nu, time, nx, ny = getFields(filename, i)
        plt.figure()
        ax = plt.gca()
        ma = np.max(omega)
        mi = np.min(omega)
        m = max(ma,-mi)
        plt.title("time = "+"{:.4f}".format(time[0]))
        im = ax.imshow(omega, extent=domain, origin="lower", cmap='seismic', vmin=-m, vmax=m, interpolation='bicubic')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.25)
        plt.colorbar(im, cax=cax, orientation="horizontal")
        plt.savefig(filename+"_omega"+str(i).zfill(5)+".png")

if __name__== "__main__":
    build_dir = "../../../build"
    app_name = os.path.join(build_dir, "apps/NavierStokes_2D_GeoRL/navier_stokes_2d_georl")
    init_name = os.path.join(build_dir, "apps/NavierStokes_2D_GeoRL/init")

    ###
    # Reynolds number is defined as
    # Re = u L / nu
    # For a Reynolds number of 100 this means nu = 0.1*.15/100
    Re=500
    ux=0.1
    radius=0.15
    nu=ux*(2*radius)/Re
    filename="test"

    itstart=1203
    iterations=itstart+20000
    final_time=100.
    for i in range(itstart,iterations):
        t = getTime(filename,i)[0];
        radiusx = radius + 0.1*np.sin(t)
        radiusy = radius - 0.1*np.sin(t)
        extrastr=""
        if i==itstart:
            extrastr=" --rv"
        command= f"{app_name} --if "+filename+"_t"+str(i)+" --of "+filename+" --T "+str(final_time)+" --origo -.5 0.025 --semimajoraxis "+str(radiusx)+" --semiminoraxis "+str(radiusy)+extrastr
        print("executing the command", command)
        tmp=execute(command).strip()
        time = getTime(filename,i+1)
        print("time=", time)

        domain, x, y, Ux, Uy, omega, nu, time, nx, ny = getFields(filename, i)

        plt.figure(1)
        plt.clf()
        ax1 = plt.gca()
        ma = np.max(omega)
        mi = np.min(omega)
        m = max(ma,-mi)
        plt.title("time = "+"{:.4f}".format(time[0]))
        im = ax1.imshow(omega, extent=domain, origin="lower", cmap='seismic', vmin=-m, vmax=m, interpolation='bicubic')
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("bottom", size="5%", pad=0.25)
        plt.colorbar(im, cax=cax, orientation="horizontal")
        plt.savefig(filename+"_omega"+str(i).zfill(5)+".png")

        #plt.figure(2)
        #plt.clf()
        #ax2 = plt.gca()
        #ma = np.max(Ux)
        #mi = np.min(Ux)
        #m = max(ma,-mi)
        #plt.title("time = "+"{:.4f}".format(time[0]))
        #im = ax2.imshow(Ux, extent=domain, origin="lower", cmap='seismic', vmin=-m, vmax=m, interpolation='bicubic')
        #divider = make_axes_locatable(ax2)
        #cax = divider.append_axes("bottom", size="5%", pad=0.25)
        #plt.colorbar(im, cax=cax, orientation="horizontal")
        #plt.savefig(filename+"_Ux"+str(i).zfill(5)+".png")

        #plt.figure(3)
        #plt.clf()
        #ax3 = plt.gca()
        #ma = np.max(Uy)
        #mi = np.min(Uy)
        #m = max(ma,-mi)
        #plt.title("time = "+"{:.4f}".format(time[0]))
        #im = ax3.imshow(Uy, extent=domain, origin="lower", cmap='seismic', vmin=-m, vmax=m, interpolation='bicubic')
        #divider = make_axes_locatable(ax3)
        #cax = divider.append_axes("bottom", size="5%", pad=0.25)
        #plt.colorbar(im, cax=cax, orientation="horizontal")
        #plt.savefig(filename+"_Uy"+str(i).zfill(5)+".png")

        if time>= final_time:
            iterations=i
            break
