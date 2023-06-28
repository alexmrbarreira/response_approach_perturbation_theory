import sys; sys.path.append('../'); 
from parameters import *
from functions import *

# This file computes the integral dq d^2P(q^2), writes it to a file and adds it to the covariance 1-loop part computed in another file

# ================================================================= #
# Define the q-Integral and compute its table
# ================================================================= #
def qIntegral(q_max):
    q_array    = klin[where(klin <= q_max)][1::]
    logq_array = log(q_array)
    integrand  = q_array**3.*Plin_int(q_array)**2.
    return integrate.trapz(integrand, logq_array)

qIntegral_table = zeros([len(k_common), len(k_common)])
for i in range(len(k_common)):
    for j in range(len(k_common)):
        k1now                 = k_common[j]
        k2now                 = k_common[i]
        qmaxnow               = min(array([fmax_r2r2*k1now, fmax_r2r2*k2now, knl_r2r2]))
        qIntegral_table[i, j] = qIntegral(qmaxnow)

covwriter_matrix(qIntegral_table, '../data_output/data_qIntegral_matrix.dat')

# ================================================================= #
# Load part of the 1-loop calculation that does not have the q-Integral and add it to them (multiply that is)
# ================================================================= #
cov_1loop_mono_noqInt = loadtxt('../data_output/data_cov_1loop_mono_noqInt.dat')
cov_1loop_quad_noqInt = loadtxt('../data_output/data_cov_1loop_quad_noqInt.dat')
cov_1loop_octu_noqInt = loadtxt('../data_output/data_cov_1loop_octu_noqInt.dat')

cov_1loop_mono = cov_1loop_mono_noqInt * qIntegral_table
cov_1loop_quad = cov_1loop_quad_noqInt * qIntegral_table
cov_1loop_octu = cov_1loop_octu_noqInt * qIntegral_table

covwriter_matrix(cov_1loop_mono, '../data_output/data_cov_1loop_mono.dat')
covwriter_matrix(cov_1loop_quad, '../data_output/data_cov_1loop_quad.dat')
covwriter_matrix(cov_1loop_octu, '../data_output/data_cov_1loop_octu.dat')
