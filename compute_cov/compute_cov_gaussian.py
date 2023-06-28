import sys; sys.path.append('../'); 
from parameters import *
from functions import *

# This file computes the angle averaged gaussian covariance:  (2/Nk)*P(k1)P(k2) delta_{12}

# Then it assumes all multipoles are the same as this angle average

# The angle dependent result is simply Cov(\vec k1, \vec k2) = P(k) if mu_12 = 1 or mu_12 = -1 (alling or anti-allign); zero otherwise

# Compute monopole gaussian covariance
cov_gauss_mono = zeros([len(k_common), len(k_common)])
for i in range(len(k_common)):
    k2_now       = k_common[i]

    # Find binwidth in k array
    index_insims = where(abs(k_sims/k2_now - 1.) == min(abs(k_sims/k2_now - 1.)))[0][0]
    if(index_insims == 0):
        binwidth_now = k_sims[index_insims+1] - k_sims[index_insims]
    elif(index_insims == len(k_sims)-1):
        binwidth_now = k_sims[index_insims] - k_sims[index_insims-1]
    else:
        binwidth_now = (k_sims[index_insims]+k_sims[index_insims+1])/2. - (k_sims[index_insims]+k_sims[index_insims-1])/2.

    index_now    = where(k_common == k2_now)[0][0]
    Vs_now       = 4.*pi*binwidth_now*k2_now**2.
    Vf           = (2.*pi)**3. / Volume
    Nk_now       = Vs_now/Vf 
    cov_gauss_mono[index_now, i] = (2./Nk_now) * Pnl_int(k2_now)**2.

covwriter_matrix(cov_gauss_mono, '../data_output/data_cov_gaussian_mono.dat')
covwriter_matrix(cov_gauss_mono, '../data_output/data_cov_gaussian_quad.dat')
covwriter_matrix(cov_gauss_mono, '../data_output/data_cov_gaussian_octu.dat')

