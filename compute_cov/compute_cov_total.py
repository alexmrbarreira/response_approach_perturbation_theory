import sys; sys.path.append('../');
from parameters import *
from functions import *

# ================================================================= #
# Load covariance multipole contributions and sum them 
# ================================================================= #
print ('squeezed_transition = ', squeezed_transition)

################### 
# Monopole
################### 

cov_1loop_mono         = loadtxt('../data_output/data_cov_1loop_mono.dat')
cov_sqtree_mono        = loadtxt('../data_output/data_cov_sqtree_mono.dat')
cov_stdtree_mono       = loadtxt('../data_output/data_cov_stdtree_mono.dat')
cov_gaussian_mono      = loadtxt('../data_output/data_cov_gaussian_mono.dat')

cov_stitched_tree_mono        = zeros([len(k_common), len(k_common)])
cov_total_mono                = zeros([len(k_common), len(k_common)])

for i in range(len(k_common)):
    k3_now       = k_common[i]
    # Build the tree level term as the sum of the standard one and the much more accurate squeezed one, in the regimes where each is valid
    cov_tree_tmp        = copy(cov_sqtree_mono[:,i])

    cov_tree_tmp[where((k_common/k3_now < squeezed_transition) & (k_common/k3_now > 1./squeezed_transition))] = cov_stdtree_mono[where((k_common/k3_now < squeezed_transition) & (k_common/k3_now > 1./squeezed_transition)), i]

    cov_stitched_tree_mono[:,i]        = cov_tree_tmp
    cov_total_mono[:,i]                = cov_gaussian_mono[:,i] + cov_stitched_tree_mono[:,i] + cov_1loop_mono[:,i]

covwriter_matrix(cov_stitched_tree_mono       , '../data_output/data_cov_stitchedtree_mono.dat')
covwriter_matrix(cov_total_mono               ,        '../data_output/data_cov_total_mono.dat')

################### 
# Quadrupole
################### 

cov_1loop_quad         = loadtxt('../data_output/data_cov_1loop_quad.dat')
cov_sqtree_quad        = loadtxt('../data_output/data_cov_sqtree_quad.dat')
cov_stdtree_quad       = loadtxt('../data_output/data_cov_stdtree_quad.dat')
cov_gaussian_quad      = loadtxt('../data_output/data_cov_gaussian_quad.dat')

cov_stitched_tree_quad        = zeros([len(k_common), len(k_common)])
cov_total_quad                = zeros([len(k_common), len(k_common)])

for i in range(len(k_common)):
    k3_now       = k_common[i]
    # Build the tree level term as the sum of the standard one and the much more accurate squeezed one, in the regimes where each is valid
    cov_tree_tmp        = copy(cov_sqtree_quad[:,i])

    cov_tree_tmp[where((k_common/k3_now < squeezed_transition) & (k_common/k3_now > 1./squeezed_transition))] = cov_stdtree_quad[where((k_common/k3_now < squeezed_transition) & (k_common/k3_now > 1./squeezed_transition)), i]

    cov_stitched_tree_quad[:,i]        = cov_tree_tmp
    cov_total_quad[:,i]                = cov_gaussian_quad[:,i] + cov_stitched_tree_quad[:,i] + cov_1loop_quad[:,i]

covwriter_matrix(cov_stitched_tree_quad       , '../data_output/data_cov_stitchedtree_quad.dat')
covwriter_matrix(cov_total_quad               ,        '../data_output/data_cov_total_quad.dat')

################### 
# Octupole
################### 

cov_1loop_octu         = loadtxt('../data_output/data_cov_1loop_octu.dat')
cov_sqtree_octu        = loadtxt('../data_output/data_cov_sqtree_octu.dat')
cov_stdtree_octu       = loadtxt('../data_output/data_cov_stdtree_octu.dat')
cov_gaussian_octu      = loadtxt('../data_output/data_cov_gaussian_octu.dat')

cov_stitched_tree_octu        = zeros([len(k_common), len(k_common)])
cov_total_octu                = zeros([len(k_common), len(k_common)])

for i in range(len(k_common)):
    k3_now       = k_common[i]
    # Build the tree level term as the sum of the standard one and the much more accurate squeezed one, in the regimes where each is valid
    cov_tree_tmp        = copy(cov_sqtree_octu[:,i])

    cov_tree_tmp[where((k_common/k3_now < squeezed_transition) & (k_common/k3_now > 1./squeezed_transition))] = cov_stdtree_octu[where((k_common/k3_now < squeezed_transition) & (k_common/k3_now > 1./squeezed_transition)), i]

    cov_stitched_tree_octu[:,i]        = cov_tree_tmp
    cov_total_octu[:,i]                = cov_gaussian_octu[:,i] + cov_stitched_tree_octu[:,i] + cov_1loop_octu[:,i]

covwriter_matrix(cov_stitched_tree_octu       , '../data_output/data_cov_stitchedtree_octu.dat')
covwriter_matrix(cov_total_octu               ,        '../data_output/data_cov_total_octu.dat')

