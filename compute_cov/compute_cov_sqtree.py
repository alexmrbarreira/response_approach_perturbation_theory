import sys; sys.path.append('../'); 
from parameters import *
from functions import *

# This file computes the covariance configuration of the tree level Trispectrum blob that is modelled as 2*\mathcal{R}_2 in the squeezed limit

#                               Cov_squeezed(k_hard)/P(k_hard)P(soft)^2 = 2 * mu_Integral[R2(k_hard, mu, -mu, -1)] / Volume

# This implementation is super poor: the angle averages are being done for every {k_i, k_j} pair, when in fact they only need to be done once. This must be improved if this is meant to be used in many generations of covariances. This was useful for some tests and then never changed back to an efficient code.

# ================================================================= #
# Define angle dependent covariance and its multipoles
# ================================================================= #

## Note: k1,k2 should be the same vector to construct a matrix that is then symmetrized.
## As is, this doesn't work well to evaluate k1=scalar, k2=another scalar

def cov_sqtree_angle(k1, k2, mu12):
    # Assume k1 is always hard and then symmetrize accordingly
    term1  = curlyR2e_nonl(k1, mu12, -mu12, -1.)
    term2  = 2.*outer(Plin_int(k2)**2., Pnl_int(k1))
    restmp = term1*term2
    return symmetrize_matrix(restmp)

def cov_sqtree_mono(k1, k2):
    # Assume k1 is always hard and then symmetrize accordingly
    term1  = (1./2)*R3e_nonl(k1) + (2./3)*R5e_nonl(k1) + (2./9)*R6e_nonl(k1) + (4./45)*R7e_nonl(k1)
    term2  = 2.*outer(Plin_int(k2)**2., Pnl_int(k1))
    restmp = term1*term2
    return symmetrize_matrix(restmp)

def cov_sqtree_quad(k1, k2):
    # Assume k1 is always hard and then symmetrize accordingly
    term1  = (2./3)*R4e_nonl(k1) + (2./9)*R6e_nonl(k1) + (8./63)*R7e_nonl(k1) 
    term2  = 2.*outer(Plin_int(k2)**2., Pnl_int(k1))
    restmp = term1*term2
    return symmetrize_matrix(restmp)

def cov_sqtree_octu(k1, k2):
    # Assume k1 is always hard and then symmetrize accordingly
    term1  = (8./35)*R7e_nonl(k1) 
    term2  = 2.*outer(Plin_int(k2)**2., Pnl_int(k1))
    restmp = term1*term2
    return symmetrize_matrix(restmp)

# ================================================================= #
# Compute moments of the squeezed covariance and write them
# ================================================================= #

cov_sqtree_mono_towrite = cov_sqtree_mono(k_common, k_common) / Volume
cov_sqtree_quad_towrite = cov_sqtree_quad(k_common, k_common) / Volume
cov_sqtree_octu_towrite = cov_sqtree_octu(k_common, k_common) / Volume

covwriter_matrix(cov_sqtree_mono_towrite, '../data_output/data_cov_sqtree_mono.dat')
covwriter_matrix(cov_sqtree_quad_towrite, '../data_output/data_cov_sqtree_quad.dat')
covwriter_matrix(cov_sqtree_octu_towrite, '../data_output/data_cov_sqtree_octu.dat')

