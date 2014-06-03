from math import pi
import atpy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from scipy.interpolate import griddata
from matplotlib import cm
#----------Declare constant values-----------------------------------
secPerYear = 3600.0 * 24.0 * 365.0
secPerDay = 24.0* 3600.0
rJup = 7.1492e9 
rSun = 6.955e10
LSun = 3.9e33
mJup = 1.89813e30
mSun = 1.99e33
AU = 1.49597871e13
gravG = 6.67259e-8

#-------Declare classes --------------------------------------
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

#------------Declare functions---------------------------------
def plot_mass_period_radius(data): # The input (data) is a pandas dataframe

    plt.rcParams['figure.figsize'] = [10, 10]  # Sets the overall size of the figure
    plt.rcParams['axes.labelsize'] = 14  # Sets the fontsize of the axis labels
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Label the axes 
    ax.set_xlabel('log10(Mass/Mjup)')
    ax.set_ylabel('Period (days)')
    ax.set_zlabel('Radius (Rjup)')

    masses = np.log10(data.Zmass.values/mJup)
    periods = map(Kepler_semi_to_period,data.ZSTA.values,data.Zmass.values,data.SEMI.values)
    radii = map(get_rinit_rfinal,data.Path.values)
    ecc_colors = map(get_ecc_color,data.ECCN.values)

    for i in range(len(masses)):
        mvals = [masses[i],masses[i]] 
        pvals = [periods[i]/secPerDay,periods[i]/secPerDay]
        rvals = radii[i]
        a = Arrow3D(mvals,pvals,rvals, mutation_scale=18, lw=1, arrowstyle="-|>", color=ecc_colors[i])
        ax.add_artist(a)
        ax.plot(mvals, pvals, rvals,'-o',color=ecc_colors[i],markersize=1)

    # Add the eccentricity colorbar to the plot
    cm = plt.cm.get_cmap('spectral')
    sc = ax.scatter([0,0], [0,0], [0,0], c=[0,200./255.], vmin=0, vmax=1, s=0, cmap=cm)
    cb = plt.colorbar(sc,shrink=0.75)
    cb.set_label('Eccentricity')
    return
#---------------------------------------------
def plot_mass_period_delta_radius(data,fignum=1): # The input (data) is a pandas dataframe
                                         # Note: all datasets need to have the SAME ECCN VALUES!!!
    masses = np.array(data.Zmass.values/mJup)
    masses = np.log10(masses)
    periods = map(Kepler_semi_to_period,data.ZSTA.values,data.Zmass.values,data.SEMI.values)
    periods = np.array(periods)/secPerDay
    periods = np.log10(periods)
    delta_radii = map(get_deltaR,data.Path.values)
    ecc_value = data.ECCN.values[0]


    #Set up and perform the interpolations...
    mass_period_input_points = np.array([masses,periods])
    mass_period_input_points = mass_period_input_points.T
    deltaR_input_values = np.array(delta_radii)
    grid_m, grid_p = np.mgrid[-2:1:20j, 0:3:20j]
    grid_deltaR = griddata(mass_period_input_points, deltaR_input_values, (grid_m, grid_p), method='nearest')
    #    print 'Successfully interpolated the data!'


    # Set up the plot
    plt.rcParams['figure.figsize'] = [5, 5]  # Sets the overall size of the figure
    plt.rcParams['axes.labelsize'] = 14  # Sets the fontsize of the axis labels
    fig = plt.figure(fignum)
    ax = fig.gca()
    ax.set_ylabel('log10 Mass (mJup)')
    ax.set_xlabel('log10 Period (days)')
    #    ax.set_xscale('log')
    #    ax.set_yscale('log')
    titlestr = 'Runs with ecc = '+str(ecc_value)
    ax.set_title(titlestr)

    # Plot the data
    p = ax.pcolor(grid_p, grid_m, grid_deltaR, cmap = cm.Spectral, vmin=-100, vmax=100)
    ax.plot(periods,masses,'ko',markersize=4)

    cb = plt.colorbar(p,shrink=0.75)
    cb.set_label('Percent change in radius')

    return 
#---------------------------------------------
def Kepler_semi_to_period(starMass,planetMass,semi):
    Mtot = starMass + planetMass
    period = ( semi**3 / (gravG * Mtot) ) **(0.5)
    period = period * 2.0 * pi

    return period # returns the period in seconds
#---------------------------------------------
def get_rinit_rfinal(filepath):
    atmos = atpy.Table(filepath+'/atmos_data.txt',type='ascii')
    rinitial = atmos.R[0]
    rfinal =atmos.R[-1]
    return [rinitial/rJup,rfinal/rJup]
#---------------------------------------------
def get_ecc_color(ecc_val):
    ecc_color = int(ecc_val * 255)
    ecc_color = plt.cm.spectral(ecc_color)
    return ecc_color
#---------------------------------------------
def get_deltaR(filepath):
    [rinit,rfinal] = get_rinit_rfinal(filepath)
    deltaR = (rfinal - rinit) / rinit
    deltaR = deltaR * 100.0
    return deltaR
#---------------------------------------------
