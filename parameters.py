from numpy import *
from scipy import *
from pylab import *
from scipy import interpolate, integrate, special, ndimage 
from numpy import linalg
import scipy.ndimage as ndimage
from matplotlib.colors import LogNorm
import pylab, os, scipy, warnings
rcParams.update({'text.usetex': False, 'mathtext.fontset': 'stix'})
warnings.filterwarnings("ignore", category=DeprecationWarning)

# ==================================================================== #
# Miscelaneous parameters
# ==================================================================== #
tinycorr  = 1.0e-10

# Put working directory here (all paths are relative to this)
localpath = '/home/barreira/a.Other/d-github_workrepos/a-response_approach_perturbation_theory/' 

# ==================================================================== #
# Needed cosmo parameters and others
# ==================================================================== #
h_hubble = 0.72 #for Blot et al 2014

Volume = 656.25**3. #This is for Blot et al 2015

# ==================================================================== #
# Power spectrum files
# ==================================================================== #

# Linear
f_pk_camb = loadtxt(localpath + 'data_input/data_powerspectra/datacamb_test_matterpower_blotetal.dat')

# Nonlinear
f_pk_coyo = loadtxt(localpath + 'data_input/data_powerspectra/data_pkemulator_blotetal.dat')

# ==================================================================== #
# Load nonlinear growth-only and full responses and create interpolators
# ==================================================================== #

f_gnl   = loadtxt(localpath + 'data_input/data_nonlinear_responses/wagner/rebin_mean_20_test10_Pk_fits_JK6_z0.dat_8') #file with G responses [k, G1, eG1, G2/2, eG2/2, etc]
k_gnl   = f_gnl[:,0]
G1_gnl  = f_gnl[:,1]
eG1_gnl = f_gnl[:,2]
G2_gnl  = f_gnl[:,3]*2.
eG2_gnl = f_gnl[:,4]*2.
G1_gnl_int = interpolate.interp1d(k_gnl, G1_gnl)
G2_gnl_int = interpolate.interp1d(k_gnl, G2_gnl)

f_rnl   = loadtxt(localpath + 'data_input/data_nonlinear_responses/wagner/rebin_mean_20_seed10_Pk_fits_JK6_z0.dat_8') #file with R responses [k, R1, eR1, R2/2, eR2/2, etc]
k_rnl   = f_rnl[:,0]
R1_rnl  = f_rnl[:,1]
eR1_rnl = f_rnl[:,2]
R2_rnl  = f_rnl[:,3]*2.
eR2_rnl = f_rnl[:,4]*2.
R3_rnl  = f_rnl[:,5]*6.
eR3_rnl = f_rnl[:,6]*6.
R1_rnl_int = interpolate.interp1d(k_rnl, R1_rnl)
R2_rnl_int = interpolate.interp1d(k_rnl, R2_rnl)
R3_rnl_int = interpolate.interp1d(k_rnl, R3_rnl)

# ==================================================================== #
# Parameters of the calculation of the 1-loop part captured by R2R2 (compute_cov/compute_cov_R2R2.py)
# ==================================================================== #
knl_r2r2  = 0.30   # k_nl h/Mpc
fmax_r2r2 = 1./2  # q is always smaller than min([fmax_r2r2*k1, fmax_r2r2*k3, knl_r2r2])

# ==================================================================== #
# Parameters of the calculation of the standard tree level covariance (compute_cov/compute_cov_stdtreelevel.py)
# ==================================================================== #
epsilon_treelevel = 1.0e-9 # this is a tiny number added to prevent mu to get to +1 or -1, which would cause numerical instabilities. Varying it shows that this is a good approximation
Nmu_treelevel     = 1.0e3 # size of mu array in angle average integrals

# ==================================================================== #
# Parameters of the calculation of the total covariance (compute_cov/compute_cov_total.py)
# ==================================================================== #
squeezed_transition = 2. # Consider accurate squeezed limit results when max{k1, k3}/min{k1, k3} > squeezed_transition


