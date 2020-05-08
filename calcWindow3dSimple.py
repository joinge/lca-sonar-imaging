"""
@file
@brief Code for computing LCA steering angles

 Sample code for computing steering angles for the
 Low Complexity Adaptive (LCA) beamformer.

 Released as complementary source code to the journal article titled
 "Low Complexity Adaptive Sonar Imaging" (Buskenes J.I, Hansen R. E., Austeng A.)
 submitted for publishing to IEEE Journal of Oceanic Engineering, 2015.

 No license restrictions applies. Use, alter and/or redistribute as you wish.
     Author:      Jo Inge Buskenes
     Affiliation: University of Oslo / The Norwegian Defense Research Establishment (FFI)
     Date:        February 2016
"""

from mynumpy import *
import mypylab as pl

  #######################################################################################
 ## WINDOW FUNCTIONS ##
######################
   
def getTrig(M=10, alpha=0.54, phi=0, d=0.5, k=2*pi):
   """ Generates a steered Trigonometric window
       
         @param[in] M         Number of array elements
         @param[in] beta      Trigonometric window parameter
         @param[in] phi       Vector of steering angles  [rad]
         @param[in] d         Element spacing            [m]
         @param[in] k         Size of wavenumber         [rad/m]

         @return              Steered Trigonometric windows
   """
   m    = linspace(0,1,M)
   wc   = alpha - (1-alpha)*cos(2*pi*m)         # Window coefficients
   m    = arange(0,M)                           # Create M indices from 0 to 1
   a    = exp(-1j*k*m*d*sin(phi))               # Steering vector
   ws   = dot(wc, a)                            # Normalization factor
   win  = a * wc / ws                           # Steered and normalized window   
   

def getKaiser(M=10, beta=3, phi=0, d=0.5, k=2*pi):
   """ Generates a steered Kaiser window
       
         @param[in] M         Number of array elements
         @param[in] beta      Kaiser window parameter
         @param[in] phi       Vector of steering angles  [rad]
         @param[in] d         Element spacing            [m]
         @param[in] k         Size of wavenumber         [rad/m]

         @return              Steered Kaiser windows
   """
   wc   = kaiser(M, beta)                       # Window coefficients
   m    = arange(0,M)                           # Create M indices from 0 to 1
   a    = exp(-1j*k*m*d*sin(phi))               # Steering vector
   ws   = dot(wc, a)                            # Normalization factor
   win  = a * wc / ws                           # Steered and normalized window
   
   return win
   
   
  #######################################################################################
 ## STEERING COMPUTATION ##
##########################
   
def find3dB(x,y):
   """ Searches for the index 'x' where 'y' is -3dB.
       Linear interpolation is used to improve the result.
       
         @param[in] x         Indexes
         @param[in] y         y=f(x), can be real of complex (absolute value used)
         @return    x_m3dB    Index where 'y' is -3 dB
   """
   N = y.shape[0]
   
   m3dB   = 1/sqrt(2)                              # -3dB
   y_abs  = abs(y)                                 # Work with the absolute value
   yr_abs = y_abs[::-1]                            # Reversed view
      
   # If the function decreases from its center onwards:
   if y_abs[N/2] > y_abs[N/2+1]:
      idx_from_center = argmax(y_abs[N/2:] < m3dB) # Search upwards for first index less than -3dB
      idx = N/2 + idx_from_center                  # Correct for the index being centered

   # If the function increases from its center onwards:
   else:
      idx_from_center = argmax(yr_abs[N/2:] < m3dB)# Search downwards for first index less than -3dB
      idx = N/2 - idx_from_center                  # Correct for the index being centered
   
   #  Improve index precision through linear interpolation
   xi, xim1 = x[idx], x[idx-1]
   yi, yim1 = y_abs[idx], y_abs[idx-1]
   idx_interpolated = xim1 + (m3dB-yim1)*(xi-xim1)/(yi-yim1)
      
   return idx_interpolated


def calc3dBSimple(M=None, lmda=None, d=None, beta=None, alpha=None, plot_file=""):
   """ Takes array and window parameters as input and returns the steering needed
       to cut the distance to the -3 dB point in half.
       
         @param[in] M           Number of array elements
         @param[in] lmda        Wavelength                      [m]
         @param[in] d           Element spacing                 [m]
         @param[in] beta        Kaiser window parameter
         @param[in] alpha       Trigonometric window parameter
         @param[in] plot_file   Filename to save the plots

         @return    phi_m3dB    Steering needed to cut distance to -3 dB point in half
   """
   if not M:    M     = 16.0                    # Number of array elements
   if not lmda: lmda  = 1.0                     # Wavelength                      [m]
   if not d:    d     = lmda / 2.0              # Element spacing                 [m]
   if beta:                                     # Kaiser window parameter.
      alpha = None                              # (takes priority over the trig. one)
   elif not alpha:
      beta  = 3
   
   D   = (M-1)*d                                # Physical length of the array    [m]
   x   = linspace(-D/2, D/2, M)                 # Element positions               [m]
   
   # Create a frequency axis that spans the 0-0 bandwidth of a rectangular window
   Nkx = 1e3                                    # Number of wavenumbers to evaluate
   kx  = linspace(-2*pi/(M*d), 2*pi/(M*d), Nkx) # Azimuth/alongtrack wavenumber   [rad/m]
   theta = arcsin(kx*lmda/(2*pi))               # Azimuth angle                   [rad]
           
   # Compute the 3dB bandwidth of a non steered window
   if   beta:  w0 = getKaiser(M=M, beta=beta,   phi=0, d=D/(M-1), k=2*pi/lmda)
   elif alpha: w0 = getTrig(  M=M, alpha=alpha, phi=0, d=D/(M-1), k=2*pi/lmda)
   W0 = dot(exp(1j*outer(kx, x)), w0)           # Window's spatial frequency response
   kx_0_m3dB    = find3dB(kx, W0)               # kx of -3dB point                [rad/m]
   theta_0_m3dB = arcsin(kx_0_m3dB*lmda/(2*pi)) # theta of -3dB point             [rad]

   # Make steering angles up 75% of the rectangular window's bandwidth
   Nphi = 1000                                   # Number of steering angles
   phi  = linspace(-arcsin(lmda/(M*d))*0.75, 0, Nphi) # Steering.                  [rad]

   # Compose a matrix of steered windows
   w = zeros((Nphi, M), dtype=complex64)
   for i in arange(phi.shape[0]):
      if   beta:  w[i] = getKaiser(M, beta=beta, phi=phi[i], d=D/(M-1), k=2*pi/lmda)
      elif alpha: w[i] = getTrig(M=M, alpha=alpha, phi=phi[i], d=D/(M-1), k=2*pi/lmda)
      
   # Compute window frequency response for each window and find its -3dB intersection
   W = dot(exp(1j*kx_0_m3dB*x/2), w.T)          # Window's spatial frequency response
   phi_m3dB = find3dB(phi, W)                   # Steering amount for -3dB point   [rad/m]

   # Plot the result (if a filename is specified)
   if plot_file != "":
      if   beta:  w_phi = getKaiser(M=M, beta=beta,   phi=phi_m3dB, d=D/(M-1), k=2*pi/lmda)
      elif alpha: w_phi = getTrig(  M=M, alpha=alpha, phi=phi_m3dB, d=D/(M-1), k=2*pi/lmda)
      W_phi = dot(exp(1j*outer(kx, x)), w_phi)
      plotW(theta, W0, W_phi, theta_0_m3dB, phi_m3dB, phi, plot_file)
   
   return phi_m3dB

# calc3dBSimple(M=16, lmda=1, d=0.5, beta=3, plot_file="test.eps")


  #######################################################################################
 ## PLOTTING ##
##############

def plotW(theta, W0, Wphi, theta_m3dB, phi_m3dB, phi, filename):
#  W0                                           # Steered response of centered window   
   theta      = theta*180/pi                    # Azimuth axis [degrees]
   phi        = phi*180/pi                      # Steering angles [degrees]
   theta_m3dB = theta_m3dB*180/pi               # Angle of -3dB point [degrees]
   phi_m3dB   = phi_m3dB*180/pi                 # Steering needed to obtain new -3dB
   theta_lims = [theta[0], theta[-1]]           # Lower and upper bound for azimuth angle
   phi_lims   = [phi[0],   phi[-1]]             # Lower and upper bound of steering
   y_lims     = [-10,5]                         # y-axis from -5 to 10 dB
   
   fn = pl.figure()
   ax = fn.add_subplot(111)
   ax.plot(theta, 20*log10(abs(W0)), 'b')       # Plot response of unsteered window
   ax.plot(theta, 20*log10(abs(Wphi)), 'r')     # Plot response of steered window

   ax.plot([theta_m3dB]*2,   y_lims, 'k')       # Mark angle where -3dB occurs
   ax.plot([theta_m3dB/2]*2, y_lims, 'k')       # Mark angle where we want the new -3dB 
   ax.plot([phi_m3dB]*2,     y_lims, 'y')       # Mark amount of steering needed
   ax.plot(theta_lims,       [20*log10(1/sqrt(2))]*2, 'k') # Mark -3dB
   ax.plot([phi_lims[0]]*2,  y_lims, 'g')       # Mark lower bound for steering angle
   ax.plot([phi_lims[1]]*2,  y_lims, 'g')       # Mark upper bound for steering angle
   ax.set_xlim(theta_lims)                      # Specify extent of x-axis
   ax.set_ylim(y_lims)                          # Specify extent of y-axis
   ax.set_xlabel(r"$\theta$ [rad]")             # Label x-axis
   ax.set_yticks([20*log10(1/sqrt(2)),0])       # Only need to mark the -3dB and 0dB
   ax.set_yticklabels(["-3","0"])               # values on the y-axis with ticks.
   ax.set_ylabel("Amplitude [dB]")              # Label y-axis
   ax.grid(True)                                # Add some grid lines
   fn.savefig(filename)                         # Save figure to file
   

  #######################################################################################
 ## USAGE EXAMPLE ##
###################

def createPhiMatrix():
   Ms = arange(64-7)+8
   ds = linspace(0.5,2.5,5)
   betas = linspace(0,5,10)
   
   phi = zeros((Ms.shape[0], ds.shape[0], betas.shape[0]))
   for m,M in enumerate(Ms):
      print M
      for di,d in enumerate(ds):
         for beta_i, beta in enumerate(betas):
            phi[m,di,beta_i] = calc3dBSimple(M=M, lmda=1, d=d, beta=beta)
   
   save('phi.npy', phi)

# createPhiMatrix()
   
