import sys; sys.path.append('../'); 
from parameters import *
from functions import *

# ================================================================= #
# Compute matrix P(ki)P(kj) to devide the covariances for plotting 
# ================================================================= #
matrix_PiPj = outer(Pnl_int(k_common), Pnl_int(k_common))
# ================================================================= #
# Load 1loop covariance with R2R2
# ================================================================= #
cov_1loop_mono = loadtxt('../data_output/data_cov_1loop_mono.dat')
cov_1loop_quad = loadtxt('../data_output/data_cov_1loop_quad.dat')
cov_1loop_octu = loadtxt('../data_output/data_cov_1loop_octu.dat')
# ================================================================= #
# Load squeezed covariance
# ================================================================= #
cov_sqtree_mono = loadtxt('../data_output/data_cov_sqtree_mono.dat')
cov_sqtree_quad = loadtxt('../data_output/data_cov_sqtree_quad.dat')
cov_sqtree_octu = loadtxt('../data_output/data_cov_sqtree_octu.dat')
# ================================================================= #
# Load standard tree level result and the stictched version
# ================================================================= #
cov_stdtree_mono = loadtxt('../data_output/data_cov_stdtree_mono.dat')
cov_stdtree_quad = loadtxt('../data_output/data_cov_stdtree_quad.dat')
cov_stdtree_octu = loadtxt('../data_output/data_cov_stdtree_octu.dat')

cov_stitchedtree_mono = loadtxt('../data_output/data_cov_stitchedtree_mono.dat')
cov_stitchedtree_quad = loadtxt('../data_output/data_cov_stitchedtree_quad.dat')
cov_stitchedtree_octu = loadtxt('../data_output/data_cov_stitchedtree_octu.dat')
# ================================================================= #
# Load standard Gaussian result
# ================================================================= #
cov_gauss_mono = loadtxt('../data_output/data_cov_gaussian_mono.dat')
cov_gauss_quad = loadtxt('../data_output/data_cov_gaussian_quad.dat')
cov_gauss_octu = loadtxt('../data_output/data_cov_gaussian_octu.dat')
# ================================================================= #
# Load total covariance, use to compute matrix C(i,i)C(j,j) 
# ================================================================= #
cov_total_mono =  loadtxt('../data_output/data_cov_total_mono.dat')
cov_total_quad =  loadtxt('../data_output/data_cov_total_quad.dat')
cov_total_octu =  loadtxt('../data_output/data_cov_total_octu.dat')

matrix_CiiCjj_mono = outer(diagonal(cov_total_mono), diagonal(cov_total_mono))
matrix_CiiCjj_quad = outer(diagonal(cov_total_quad), diagonal(cov_total_quad))
matrix_CiiCjj_octu = outer(diagonal(cov_total_octu), diagonal(cov_total_octu))
