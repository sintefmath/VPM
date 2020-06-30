import numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import cm
import math
from matplotlib import rc

from scipy import signal

from scipy.fftpack import fft

#prec = np.float32
prec = np.float64

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

def readConf(filename, timestep, vel=""):
    if True:
        conffilename = filename+"_t"+str(timestep)+"_conf"+vel+".dat"
        llx, lly, urx, ury, eps, nu, time, sigma, population_threshold, bc_treshold_xl, bc_treshold_xr, bc_treshold_yl, bc_treshold_yr, Uinfty_x, Uinfty_y = readFromFile(conffilename, None, prec, 15)
        seek_num = 15
    else:# old data files without Uinfty
        conffilename = filename+"_t"+str(timestep)+"_conf"+vel+".dat"
        llx, lly, urx, ury, eps, nu, time, sigma, population_threshold, bc_treshold_xl, bc_treshold_xr, bc_treshold_yl, bc_treshold_yr = readFromFile(conffilename, None, prec, 13)
        seek_num = 13

    if prec==np.float32:
        seek_num *= 4
    if prec==np.float64:
        seek_num *= 8

    remesh_ison, remesh_steps , remesh_method, N, nx, ny, order_ODEsolver, bc_xl, bc_xr, bc_yl, bc_yr = readFromFile(conffilename, None, np.int32, -1, seek_num)

    domain = np.array((llx, urx, lly, ury))
    remesh_ison = bool(remesh_ison)
    remesh_method = int(remesh_method)

    return domain, nx, ny, N, nu, time



def readParticles(filename, timestep, velocityField=False, backgroundField=False):

    vel_str = ""
    #if velocityField:
    #    vel_str = "_vel"

    domain, nx, ny, N, nu, time = readConf(filename, timestep, vel_str)

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

def plotSideWind():
    ax = plt.axes()
    L=1
    U=.25
    LL=L*L
    for x in np.linspace(-.95,1.,25,endpoint=False):
        pos=x+1.
        if pos<L:
            wind=2.*U/LL*(L*pos - pos*pos/2.)
        else:
            wind=U
        ax.arrow(-1., x, wind, 0, head_width=0.05, head_length=0.1, fc='k', ec='k')

def relErrorTime(filename, timesteps):
    L1_error = []
    L2_error = []
    Linfty_error = []
    times = []
    for t in np.arange(timesteps):
        domain, x, y, omega, nu, time, nx, ny = readParticles(filename, t)
        r2 = x**2 + y**2
        strength=1.
        rct_squared = 0.15**2+4.*nu*time;
        omega0 = strength/(math.pi*rct_squared)*np.exp(-r2/rct_squared)

        #L2_error.append(np.sqrt(np.sum(np.multiply(omega-omega0, omega-omega0)))/np.sqrt(np.sum(np.multiply(omega0,omega0))))
        #Linfty_error.append(np.max(np.abs(omega-omega0))/np.max(np.abs(omega0)))
        L1_error.append(np.sum(np.abs(omega-omega0))/np.sum(np.abs(omega0)))
        times.append(time)
    domain, nx, ny, N, nu, time = readConf(filename, 0)
    plt.ion()
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.title(r'Relative error for the VPM with $\nu=$%s'%nu)
    plt.plot(times,L1_error,'-x', label=r'$L_1$-error grid=%sx%s'%(nx,ny))
    #plt.plot(times,L2_error,'-x', label=r'$L_2$-error grid=%sx%s'%(nx,ny))
    #plt.plot(times,Linfty_error,'--x', label=r'$L_\infty$-error grid=%sx%s'%(nx,ny))
    plt.xlabel('time')
    plt.ylabel('error')
    #plt.legend(loc='best')
    plt.show()

def relError(filename, timestep, target_time=1.):
    time=-1
    if timestep==-1:
        timestep=0
        found=False
        while abs(time-target_time)>1e-6 :
            domain, x, y, omega, nu, time, nx, ny = readParticles(filename, timestep)
            timestep += 1
        if abs(time-target_time)>1e-6 :
            print("Time step could not be found")
            return
    else:
        domain, x, y, omega, nu, time, nx, ny = readParticles(filename, timestep)

    #print(nu, time)
    r2 = x**2 + y**2
    strength=1.
    rct_squared = 0.15**2+4.*nu*time;
    omega0 = strength/(math.pi*rct_squared)*np.exp(-r2/rct_squared)

    #print(np.max(np.max(omega0)))

    L1_error = np.sum(np.abs(omega-omega0))/np.sum(np.abs(omega0))
    L2_error = np.sqrt(np.sum(np.multiply(omega-omega0, omega-omega0)))/np.sqrt(np.sum(np.multiply(omega0,omega0)))
    Linfty_error = np.max(np.abs(omega-omega0))

    #print(np.max(np.abs(omega)), np.max(np.abs(omega0)))

    #plt.ion()
    matplotlib.rcParams.update({'font.size': 22})

    plt.figure(1)
    plt.clf()
    plt.subplot(311)
    plt.title("Error t=%s"%time)
    plt.imshow(np.abs(omega-omega0), extent=domain, origin="lower")#, cmap=blue_red1)
    plt.colorbar()

    plt.subplot(312)
    plt.title("VPM t=%s"%time)
    plt.imshow(omega, extent=domain, origin="lower")#, cmap=blue_red1)
    plt.clim([0,14.022])
    plt.colorbar()

    plt.subplot(313)
    plt.title("Exact t=%s"%time)
    plt.imshow(omega0, extent=domain, origin="lower")#, cmap=blue_red1)
    plt.clim([0,14.022])
    plt.colorbar()
    

    plt.figure(3)
    plt.clf()
    plt.title("Error t=%s"%time)
    plt.imshow(np.abs(omega-omega0), extent=domain, origin="lower")#, cmap=blue_red1)
    plt.clim([0,0.01])
    plt.colorbar()
    plt.savefig(filename+"_t"+str(timestep).zfill(4)+"_relerror_omega_m_exact.png")

    plt.figure(3)
    plt.clf()
    plt.title("Exact t=%s"%time)
    plt.imshow(omega0, extent=domain, origin="lower")#, cmap=blue_red1)
    plt.clim([0,14.022])
    plt.colorbar()
    plt.savefig(filename+"_t"+str(timestep).zfill(4)+"_relerror_exact.png")

    plt.figure(4)
    plt.clf()
    plt.title("VPM t=%s"%time)
    plt.imshow(omega, extent=domain, origin="lower")#, cmap=blue_red1)
    plt.clim([0,14.022])
    plt.colorbar()
    plt.savefig(filename+"_t"+str(timestep).zfill(4)+"_relerror_omega.png")

    return time, timestep, L1_error, L2_error, Linfty_error

def plotParticles(filename='', timestep=0, regular = True, save = False, case="nothing"):
    masked = False

    if not(save):
        plt.ion()


    if False:
        domain, x, y, Ux, Uy, nu, time, nx, ny = readParticles(filename, timestep, True)
        im = plt.scatter(x, y, c=Ux, lw=0)#, cmap=blue_red1)
    else:
        domain, x, y, omega, nu, time, nx, ny = readParticles(filename, timestep)

        if regular:
            #plt.imshow(omega, interpolation='none', extent=domain, origin="lower")
            if masked:
                omega_masked = np.ma.masked_where(abs(omega) < 1e-6, omega)
                im = plt.imshow(omega_masked, extent=domain, origin="lower")#, vmin=-40, vmax=40)#, cmap=blue_red1)
            else:
                im = plt.imshow(omega, extent=domain, origin="lower")#, vmin=-40, vmax=40)#, cmap=blue_red1)
        else:
            im = plt.scatter(x, y, c=omega, lw=0)#, cmap=blue_red1)

    plt.gca().set_aspect('equal', adjustable='box')

    #plt.xlim(-1.0, 2.0)
    #plt.ylim(-2.0, 1.0)

    plt.title(r'time $t=$%.2f'%time)

    XX = np.linspace(domain[0], domain[1], nx)

    if case=="hill":
### hill
        plt.plot(XX,-.75+.1*np.exp(-.1*((XX-.25)/0.05)**2),'k')
        cl = min(abs(np.max(omega)),abs(np.min(omega)))
        plt.clim([-cl,cl])
    elif case=="circle":
### circle
        theta = np.linspace(-np.pi, np.pi, 200)
        plt.plot(np.sin(theta), np.cos(theta),'k')
        plt.xlim([-1.1,1.1])
        plt.ylim([-1.1,1.1])
        #plt.clim([-10,10])
    elif case=="bluff":
### circle bluff body
        CIRCLE_R = .2#-(x[0,1]-x[0,0])
        CIRCLE_X = -.5
        CIRCLE_Y = .111
        theta = np.linspace(-np.pi, np.pi, 200)
        plt.plot(CIRCLE_X+CIRCLE_R*np.sin(theta), CIRCLE_Y+CIRCLE_R*np.cos(theta),'k')
        #plt.xlim([-.8,2.2])
        #plt.ylim([-.5,.5])
        plt.xlim([-.8,3.2])
        plt.ylim([-1.,1.])
        #plt.clim([-6,6])
        #plt.clim([-2,2])
    elif case=="bluffhalf":
### circle bluff body
        CIRCLE_R = .2#-(x[0,1]-x[0,0])
        CIRCLE_X = -.5
        CIRCLE_Y = .111
        theta = np.linspace(-np.pi, 0, 200)
        plt.plot(CIRCLE_X+CIRCLE_R*np.sin(theta), CIRCLE_Y+CIRCLE_R*np.cos(theta),'k')
        plt.plot(np.array([CIRCLE_X, CIRCLE_X]), np.array([CIRCLE_Y-CIRCLE_R, CIRCLE_Y+CIRCLE_R]),'k')
        #plt.xlim([-1.,2.])
        #plt.ylim([-.75,.75])
    elif case=="flat":
### flat
        plt.plot(XX,-.75+.0*XX,'k')
    elif case=="ellipse":
### ellipse
        ELLIPSE_R = .2#-(x[0,1]-x[0,0])
        ELLIPSE_X = -.4
        ELLIPSE_Y = .011
        DEGREE=45.
        angle = 2.*np.pi*DEGREE/360.
        A=.25
        B=1.25
        #
        theta = np.linspace(-np.pi, np.pi, 200)
        xd = + A*np.cos(theta)*np.cos(angle) - B*np.sin(theta)*np.sin(angle)
        yd = + A*np.cos(theta)*np.sin(angle) + B*np.sin(theta)*np.cos(angle)
        plt.plot(ELLIPSE_X+ELLIPSE_R*xd, ELLIPSE_Y+ELLIPSE_R*yd,'k')
        plt.xlim([-.8,3.2])
        plt.ylim([-1.,1.])
        #plt.clim([-6,6])
        #plt.clim([-2,2])
        plt.clim([-1,1])



    #scale colorbar to the height of the image
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    plt.tight_layout(pad=0, w_pad=0, h_pad=0)
    if save:
        plt.savefig(filename+"_t"+str(timestep).zfill(4)+"_omega.png")
    else:
        plt.show()

def plotVelocity(filename='', timestep=0, streamplot = False, backgroundField = False, save = False):

    if backgroundField:
        domain, x, y, Ux, Uy, Ux_bg, Uy_bg, nu, time, nx, ny = readParticles(filename, timestep, True, True)
    else:
        domain, x, y, Ux, Uy, nu, time, nx, ny = readParticles(filename, timestep, True)

    if not(save):
        plt.ion()

    if Ux.ndim==1:
        plt.figure(1)
        plt.clf()
        plt.scatter(x, y, c=Ux, lw=0)#, cmap=blue_red1)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.colorbar()

        plt.figure(2)
        plt.clf()
        plt.scatter(x, y, c=Uy, lw=0)#, cmap=blue_red1)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.colorbar()

        if backgroundField:
            plt.figure(3)
            plt.clf()
            plt.scatter(x, y, c=Ux_bg, lw=0)#, cmap=blue_red1)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.colorbar()

            plt.figure(4)
            plt.clf()
            plt.scatter(x, y, c=Uy_bg, lw=0)#, cmap=blue_red1)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.colorbar()
    else:
        speed = np.sqrt(Ux**2 + Uy**2)
        if backgroundField:
            speed_bg = np.sqrt(Ux_bg**2 + Uy_bg**2)
        if streamplot:
            filename = filename+"_stream"
            plt.streamplot(x,y,Ux, Uy,
                        color=speed,         # array that determines the colour
                        cmap=cm.cool,        # colour map
                        linewidth=2,         # line thickness
                        arrowstyle='->',     # arrow style
                        arrowsize=1.5)       # arrow size
            plt.xlim(-1.0, 1.0)
            plt.ylim(-1.0, 1.0)
            plt.colorbar()                      # add colour bar on the right
            if backgroundField:
                filename = filename+"_stream"
                plt.streamplot(x,y,Ux_bg, Uy_bg,
                            color=speed_bg,         # array that determines the colour
                            cmap=cm.cool,        # colour map
                            linewidth=2,         # line thickness
                            arrowstyle='->',     # arrow style
                            arrowsize=1.5)       # arrow size
                plt.xlim(-1.0, 1.0)
                plt.ylim(-1.0, 1.0)
                plt.colorbar()                      # add colour bar on the right
        else:
            filename = filename+"_quiver"
            plt.quiver(x,y,Ux, Uy)
            plt.xlim(-1.0, 1.0)
            plt.ylim(-1.0, 1.0)

        plt.figure(4)
        plt.imshow(speed, extent=domain, origin="lower")
        plt.colorbar()                      # add colour bar on the right
        plt.xlim(-1.0, 1.0)
        plt.ylim(-1.0, 1.0)

        plt.figure(5)
        plt.imshow(Ux, extent=domain, origin="lower")
        plt.colorbar()                      # add colour bar on the right
        plt.xlim(-1.0, 1.0)
        plt.ylim(-1.0, 1.0)

        plt.figure(6)
        plt.imshow(Uy, extent=domain, origin="lower")
        plt.colorbar()                      # add colour bar on the right
        plt.xlim(-1.0, 1.0)
        plt.ylim(-1.0, 1.0)
        if backgroundField:
            plt.figure(7)
            plt.imshow(speed_bg, extent=domain, origin="lower")
            plt.colorbar()                      # add colour bar on the right
            plt.xlim(-1.0, 1.0)
            plt.ylim(-1.0, 1.0)

            plt.figure(8)
            plt.imshow(Ux_bg, extent=domain, origin="lower")
            plt.colorbar()                      # add colour bar on the right
            plt.xlim(-1.0, 1.0)
            plt.ylim(-1.0, 1.0)

            plt.figure(9)
            plt.imshow(Uy_bg, extent=domain, origin="lower")
            plt.colorbar()                      # add colour bar on the right
            plt.xlim(-1.0, 1.0)
            plt.ylim(-1.0, 1.0)

        if save:
            plt.savefig(filename+"_t"+str(timestep).zfill(4)+"_U.png")
        else:
            plt.show()


def plotParticles_all(filename='', timesteps=10):
    for a in range(timesteps):
        plt.figure(a)
        particles2d(filename,a)
        plt.xlim(-1.5, 1.5)
        plt.ylim(-1.5, 1.5)

def saveAllParticles(filename, ts, te, regular, case="nothing"):
    for t in range(ts, te):
        print(t)
        plt.clf()
        plotParticles(filename, t, regular, True, case)

def timeSeriesOmega(filename, ts, te, px, py):
    timeseries_t = np.zeros(te-ts)
    timeseries_omega = np.zeros(te-ts)
    c=-1
    for t in range(ts, te):
        c += 1
        #print(t)
        domain, x, y, omega, nu, time, nx, ny = readParticles(filename, t)
        timeseries_t[c] = time
        print(x[py,px],y[py,px])
        timeseries_omega[c] = omega[py,px]
    np.save("timeseries_t_"+filename,timeseries_t)
    np.save("timeseries_omega_"+filename,timeseries_omega)
    return timeseries_t, timeseries_omega

def timeSeriesUx(filename, ts, te, px, py):
    timeseries_t = np.zeros(te-ts)
    timeseries_u = np.zeros(te-ts)
    c=-1

    for t in range(ts, te):
        c += 1
        domain, x, y, Ux, Uy, nu, time, nx, ny = readParticles(filename, t, True)
        timeseries_t[c] = time
        found = False
        dx=(domain[1]-domain[0])/nx
        for i in np.arange(Ux.shape[0]):
            rx = (x[i]-px)
            ry = (y[i]-py)
            if rx**2 + ry**2 < 0.5*dx**2:
                timeseries_u[c] = Ux[i]
                found = True
                break
        if not found:
            print(c, "not found!")
            timeseries_u[c] = -10
    np.save("timeseries_t_"+filename,timeseries_t)
    np.save("timeseries_u_"+filename,timeseries_u)
    plt.ion()
    plt.clf()
    plt.plot(timeseries_t, timeseries_u)
    plt.savefig("timeseries_"+filename+".png")
    return timeseries_t, timeseries_u


def fourier(t,w,ts,te,n):
    from scipy.interpolate import interp1d
    f2 = interp1d(t, w, kind='cubic')

    x = np.linspace(ts, te, num=n, endpoint=True)
    y = f2(x)

    xdft = np.abs(np.fft.fft(y,n))**2
    xdftabs = xdft
    freq = np.arange(n)
    maxAmp = np.max(xdftabs)
    maxAmpFreq = freq[xdftabs>freq-0.0001]
    print("maxAmp", maxAmp, maxAmpFreq)

    plt.ion()
    plt.figure(1)
    plt.plot(t,w,'kx-')
    plt.plot(x,y,'ro-')
    #plt.plot(t,0.5*np.sin(maxAmpFreq[1]*math.pi*t/10.),'b')
    #plt.plot(t,0.5*np.sin(maxAmpFreq[2]*math.pi*t/10.),'g')

    plt.plot(x,0.5*np.sin(maxAmpFreq[1]*math.pi*x/10.),'b')
    plt.plot(x,0.5*np.sin(maxAmpFreq[2]*math.pi*x/10.),'g')

    plt.figure(2)
    plt.plot(freq, xdft)


def spectrumnot(t,w,ts,te,n):
    from scipy.interpolate import interp1d
    f2 = interp1d(t, w, kind='cubic')

    x = np.linspace(ts, te, num=n, endpoint=True)
    y = f2(x)

    #from __future__ import division

    ps = np.abs(np.fft.fft(y))**2

    time_step = x[1]-x[0]
    freqs = np.fft.fftfreq(ps.size, time_step)
    idx = np.argsort(freqs)

    plt.ion()
    plt.plot(freqs[idx], ps[idx])

def tee():
    n=100
    t = np.linspace(0,1,n, endpoint=True)
    x = 1*np.cos(2.*math.pi*t)# + 1.*np.sin(8*2*math.pi*t)
    #xdft = np.fft.fft(x,n)
    #xdftabs = np.abs(xdft)
    #freq = np.arange(n)
    #maxAmp = np.max(xdftabs)
    #maxAmpFreq = freq[xdftabs>freq-0.01]
    #print("maxAmp", maxAmp, maxAmpFreq)

    N=2048
    x = np.abs(np.fft.fft(x,N))
    x = np.fft.fftshift(x)

    F = np.linspace(-N/2,N/2-1,N)/N

    plt.ion()
    plt.plot(F,x)

    #freqRel = freq[1:round(n/2)]
    #xdftRel = np.abs(xdft[1:round(n/2)])

    #plt.ion()
    #plt.plot(freqRel, xdftRel,'k')
    #print(freqRel.shape, xdftRel.shape)
    #plt.figure()
    #plt.plot(t,x)


def spectrum(t,w,ts,te,n):
    #Fs = 1000 # sampling frequency 1 kHz
    #t = 0 : 1/Fs : 0.296; % time scale
    #f = 200; % Hz, embedded dominant frequency
    #x = cos(2*pi*f*t) + randn(size(t)); % time series
    #plt.plot(t,x), axis(’tight’), grid(’on’), title(’Time series’), figure
    #nfft = 512; % next larger power of 2

    from scipy.interpolate import interp1d
    f2 = interp1d(t, w, kind='cubic')

    t = np.linspace(ts, te, num=n, endpoint=True)
    #x = np.sin(2*math.pi*t)
    x = f2(t)
    x = signal.detrend(x)#remove linear trends

    yf = fft(x) # Fast Fourier Transform

    T = t[1]-t[0]
    N=len(t)
    xf = np.linspace(0.0, 1.0/(2.0*T), int(N/2))
    plt.ion()
    plt.clf()
    plt.plot(xf, 2.0/N * np.abs(yf[0:int(N/2)]**2))
    plt.xlim(.0, .25)
    plt.ylim(.0, .04)

    plt.savefig("powspectrum_now.png")

    #y = np.abs(y**2) # raw power spectrum density

   #' plt.ion()
   #' plt.figure(1)
   #' plt.clf()
   #' plt.plot(t,x)

   #' #Y = np.fft.fft(X)

   #' L=y.shape[0]
   #' Fs=L

   #' P2 = np.abs(y/L)
   #' P1 = P2[1:L/2+1]
   #' P1[2:-2] = 2*P1[2:-2]
   #' f = Fs*np.arange(0,int(L/2))/L
   #' plt.figure(2)
   #' plt.clf()
   #' plt.plot(f,P1)

    #plt.figure(2)
    #plt.clf()
    #plt.plot(y)
    #nfft = y.shape[0]
    #h = 1+int(nfft/2)
    ##y = y[:h] # half-spectrum
    ##print("maximum frequency=", np.max(y)) # find maximum
    ##plt.ion()
    ##plt.plot(np.arange(h),y)
    ##plt.figure()
    ##plt.plot(t,x)
    #f scale = (0:nfft/2)* Fs/nfft; # frequency scale
    #plot(f scale, y),axis(’tight’),grid(’on’),title(’Dominant Frequency’)
    #f est = f scale(k); # dominant frequency estimate
    #fprintf(’Dominant freq.: true %f Hz, estimated %f Hz\n’, f, f est)
    #fprintf(’Frequency step (resolution) = %f Hz\n’, f scale(2))


def final():
    n = np.arange(150)/9.
    x1 = np.cos(2*math.pi*n/10.)
    N= 2048
    x = np.abs(np.fft.fft(x1,N))
    x = np.fft.fftshift(x)

    F = np.linspace(-N/2,N/2-1,N)/N

    plt.ion()
    plt.plot(F,x)

    plt.figure()
    plt.plot(n,x1)
def finals():
    n = np.arange(150)
    x1 = np.cos(2*math.pi*n/10.)
    N= 2048
    x = np.abs(np.fft.fft(x1,N))
    x = np.fft.fftshift(x)

    F = np.linspace(-N/2,N/2-1,N)/N

    plt.ion()
    plt.plot(F,x)

    plt.figure()
    plt.plot(n,x1)

def ddd(Fs):
    #Fs = 1000            # Sampling frequency
    T = 1./Fs             # Sampling period
    #L = 1000             # Length of signal
    L = Fs
    t = np.arange(0,L-1)*T        # Time vector

    S = 0.7*np.sin(2*math.pi*50*t) + np.sin(2*math.pi*120*t)

    X = S + 2*np.random.randn(t.shape[0])

    plt.ion()
    plt.figure(1)
    plt.clf()
    plt.plot(t,X)

    Y = np.fft.fft(X)

    P2 = np.abs(Y/L)
    P1 = P2[1:L/2+1]
    P1[2:-2] = 2*P1[2:-2]

    f = Fs*np.arange(0,int(L/2))/L
    plt.figure(2)
    plt.clf()
    plt.plot(f,P1)
