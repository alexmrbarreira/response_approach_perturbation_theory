import sys; sys.path.append('../');
from parameters import *
from functions import *

# This file computes the 1loop part of the covariance that is modelled as the product of two second-order responses

#                               Cov_ng_1loop(k1,k2)/q-Integral q^2P(q)^2 = P(k1)P(k2)/(2pi)^2 * 
#                                                                          mu_Integral R2(k1, mu, -mu, -1)*R2(k2, mu, -mu, -1)
# The q-integral is added in a separate file for speed

# ================================================================= #
# Define the mu-Integral (see notes for meaning of terms)
# ================================================================= #
# Auxiliary term definitions
def Aterm(k): return (1./2)*R3e_nonl(k) + (2./3)*R5e_nonl(k) + (2./9)*R6e_nonl(k) 
def Bterm(k): return (2./3)*R4e_nonl(k) + (2./9)*R6e_nonl(k)
def Cterm(k): return (4./9)*R7e_nonl(k)

def muIntegral_l(k1, k2, mu12):
    # Define k-dependent functions
    A1 = Aterm(k1)
    B1 = Bterm(k1)
    C1 = Cterm(k1)
    A2 = Aterm(k2)
    B2 = Bterm(k2)
    C2 = Cterm(k2)
    # Define a few terms
    term1 = 2.*A1*A2
    term2 = (2./5)*B1*B2*Plegendre_2(mu12)
    term3 = (2./5)*(A1*C2 + A2*C1)
    term4 = (4./35)*(B1*C2 + B2*C1)*Plegendre_2(mu12)
    term5 = (2./35)*(1. + 2.*Plegendre_2(mu12)**2.)*C1*C2
    return term1 + term2 + term3 + term4 + term5

# ================================================================= #
# Define the Covariance(k1,k2,mu12) function multipoles
# Note, they do not include the q-Integral term. This is done separately in another file, which helps to speed up variations of q_max
# ================================================================= #
def cov_1loop_angle(k1, k2, mu12):
    term1 = Pnl_int(k1)*Pnl_int(k2)/(2.*pi)**2.
    term3 = muIntegral_l(k1, k2, mu12)
    return 2. * term1 * term3

def cov_1loop_mono_noqInt(k1, k2):
    # Define k-dependent functions
    A1 = Aterm(k1)
    B1 = Bterm(k1)
    C1 = Cterm(k1)
    A2 = Aterm(k2)
    B2 = Bterm(k2)
    C2 = Cterm(k2)
    term1 = Pnl_int(k1)*Pnl_int(k2)/(2.*pi)**2.
    term3 = 2.*A1*A2 + (2./5)*(A1*C2 + A2*C1) + (2./25)*C1*C2
    return 2. * term1 * term3

def cov_1loop_quad_noqInt(k1, k2):
    # Define k-dependent functions
    A1 = Aterm(k1)
    B1 = Bterm(k1)
    C1 = Cterm(k1)
    A2 = Aterm(k2)
    B2 = Bterm(k2)
    C2 = Cterm(k2)
    term1 = Pnl_int(k1)*Pnl_int(k2)/(2.*pi)**2.
    term3 = (2./5)*B1*B2 + (4./35)*(B1*C2 + B2*C1) + (8./245)*C1*C2
    return 2. * term1 * term3

def cov_1loop_octu_noqInt(k1, k2):
    # Define k-dependent functions
    A1 = Aterm(k1)
    B1 = Bterm(k1)
    C1 = Cterm(k1)
    A2 = Aterm(k2)
    B2 = Bterm(k2)
    C2 = Cterm(k2)
    term1 = Pnl_int(k1)*Pnl_int(k2)/(2.*pi)**2.
    term3 = (72./1225)*C1*C2
    return 2. * term1 * term3

# ================================================================= #
# Compute and write covariance matrix multipoles
# ================================================================= #
cov_1loop_mono_noqInt_towrite = zeros([len(k_common), len(k_common)])
cov_1loop_quad_noqInt_towrite = zeros([len(k_common), len(k_common)])
cov_1loop_octu_noqInt_towrite = zeros([len(k_common), len(k_common)])
for i in range(len(k_common)):
    if(mod(i, int(len(k_common)/10))==0):
        print ('Up to line', i, 'out of', len(k_common),' ::  Each part gets faster though ... ')
    k2now = k_common[i]
    for j in range(i, len(k_common)):
        k1now          = k_common[j]
        cov_1loop_mono_noqInt_towrite[i,j] = cov_1loop_mono_noqInt(k1now, k2now)
        cov_1loop_quad_noqInt_towrite[i,j] = cov_1loop_quad_noqInt(k1now, k2now)
        cov_1loop_octu_noqInt_towrite[i,j] = cov_1loop_octu_noqInt(k1now, k2now)

# Symmetrize for k1 <--> k2
cov_1loop_mono_noqInt_towrite = symmetrize_matrix(cov_1loop_mono_noqInt_towrite) / Volume
cov_1loop_quad_noqInt_towrite = symmetrize_matrix(cov_1loop_quad_noqInt_towrite) / Volume
cov_1loop_octu_noqInt_towrite = symmetrize_matrix(cov_1loop_octu_noqInt_towrite) / Volume

covwriter_matrix(cov_1loop_mono_noqInt_towrite, '../data_output/data_cov_1loop_mono_noqInt.dat')
covwriter_matrix(cov_1loop_quad_noqInt_towrite, '../data_output/data_cov_1loop_quad_noqInt.dat')
covwriter_matrix(cov_1loop_octu_noqInt_towrite, '../data_output/data_cov_1loop_octu_noqInt.dat')

