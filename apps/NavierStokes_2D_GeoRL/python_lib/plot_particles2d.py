import numpy as np
import time
import pyVPM
import tqdm
# from mpi4py import MPI

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
from os import path, read

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

def getFields(pf):

    domain = pf.params.domain()
    nx = pf.params.nx()
    ny = pf.params.ny()
    nu = pf.params.nu()
  
    reshape = (ny, nx)
    prec=np.float64
     
    x = np.reshape(np.array([point.x for point in pf.positions]), reshape)
    y = np.reshape(np.array([point.y for point in pf.positions]), reshape)
    
    Ux = np.reshape(np.array([vel.x for vel in pf.velocity]), reshape)
    Uy = np.reshape(np.array([vel.y for vel in pf.velocity]), reshape)
    omega = np.reshape(np.array(pf.omega), reshape)
    time = pf.time
    return domain, x, y, Ux, Uy, omega, nu, time, nx, ny


def make_particle_field(nu, domain, Uinfty, nx=100, ny=50, order=1):

    dt = 0.01
  


    remesh_isOn = True
    remesh_steps = 1
    remesh_method = "M6d_fast"


    population_threshold = -1

    bc_xl = "none"
    bc_to_xl = 0.0
    bc_xr = "none"
    bc_to_xr = 0.0
    bc_yl = "none"
    bc_to_yl = 0.0
    bc_yr = "none"
    bc_to_yr = 0.0

    remeshParams = pyVPM.RemeshParams(
            remesh_isOn,
            remesh_steps,
            remesh_method
    )

    # std::shared_ptr<VPM::BCParams> bcParams = std::make_shared<VPM::BCParams>(
    #         bc_xl,
    #         bc_to_xl,
    #         bc_xr,
    #         bc_to_xr,
    #         bc_yl,
    #         bc_to_yl,
    #         bc_yr,
    #         bc_to_yr
    #         );

    # VPM::Parameters params = VPM::Parameters(
    #         domain_ll,
    #         domain_ur,
    #         m_numberofParticlesx,
    #         m_numberofParticlesy,
    #         nu,
    #         population_threshold,
    #         remeshParams,
    #         bcParams,
    #         order,
    #         Uinfty
    #         );

def run_one_step(vpm_instance, T, origo, semimajoraxis, semiminoraxis):
    random_velocity_dist = False


if __name__== "__main__":
    ###
    # Reynolds number is defined as
    # Re = u L / nu
    # For a Reynolds number of 100 this means nu = 0.1*.15/100
    Re=500
    ux=0.1
    radius=0.15
    nu=ux*(2*radius)/Re
    filename="test"
    pyVPM.initialize_module(sys.argv)
    vpm_instance = pyVPM.VPM2d(sys.argv)

    
    pf = pyVPM.ParticleField()
    itstart=1203
    iterations=itstart+20000
    final_time=100.
    random_velocity_dist = True
    center = pyVPM.Point2d(-0.5, 0.025)
    structure = pyVPM.Structure_Ellipse(center, radius, radius)
    vpm_instance.setStructure(structure)

    pyVPM.readParticlesFromFile(filename + "_t" + str(itstart), pf, random_velocity_dist)
   
    time_spent_in_solving = 0
    start_loop = time.time()
    for i in tqdm.tqdm(range(itstart,iterations)):
        # print(f"{i=} ({(i-itstart)/(iterations-itstart)*100} %)")
        t = pf.time
        radiusx = radius + 0.1*np.sin(t)
        radiusy = radius - 0.1*np.sin(t)
            
        structure = pyVPM.Structure_Ellipse(center, radiusx, radiusy)
        vpm_instance.setStructure(structure)
        start_time = time.time()
        vpm_instance.run_without_writer(pf, final_time, True)
        end_time = time.time()
        time_spent_in_solving += end_time - start_time
        

        domain, x, y, Ux, Uy, omega, nu, t, nx, ny = getFields(pf)

        # plt.figure(1)
        # plt.clf()
        # ax1 = plt.gca()
        # ma = np.max(omega)
        # mi = np.min(omega)
        # m = max(ma,-mi)
        # plt.title("time = "+"{:.4f}".format(time))
        # im = ax1.imshow(omega, extent=domain, origin="lower", cmap='seismic', vmin=-m, vmax=m, interpolation='bicubic')
        # divider = make_axes_locatable(ax1)
        # cax = divider.append_axes("bottom", size="5%", pad=0.25)
        # plt.colorbar(im, cax=cax, orientation="horizontal")
        # plt.savefig(filename+"_omega"+str(i).zfill(5)+".png")

        
        if t>= final_time:
            iterations=i
            break
    end_loop = time.time()

    print(f"Time spent in loop: {end_loop - start_loop} s")
    print(f"Time spent in solving: {time_spent_in_solving} s")
