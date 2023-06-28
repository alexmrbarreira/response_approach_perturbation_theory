import sys; sys.path.append('../'); from parameters import *
import sys; sys.path.append('../'); from functions import *

# This file computes the standard tree level covariance:  (1/2) int d\mu T(k1, -k1, k3, -k3)

# This calculation does not make the average along the radial direction of the spherical k-shells. Instead, it assumes that the result is independent of the bin width, and implicitly takes it to be sufficiently small that one can ignore the integral

# ================================================================= #
# Define the angle-dependent covariance 
# ================================================================= #

# Define the integrand
def cov_stdtree_angle(k1, k2, mu12):
    # Shorthand notation of vector diferences and angles
    k1mk2          = sqrt(k1**2. + k2**2. - 2.*k1*k2*mu12) # |\vec{k1}-\vec{k2}|
    k1pk2          = sqrt(k1**2. + k2**2. + 2.*k1*k2*mu12) # |\vec{k1}+\vec{k2}|
    ang_k1mk2andk2 = (k1*mu12-k2)/k1mk2                    # Cosine angle between \vec{k1}-\vec{k2} and \vec{k2}
    ang_k1pk2andk2 = (k1*mu12+k2)/k1pk2                    # Cosine angle between \vec{k1}+\vec{k2} and \vec{k2}
    ang_k2mk1andk1 = (k2*mu12-k1)/k1mk2                    # Cosine angle between \vec{k2}-\vec{k1} and \vec{k1}
    ang_k2pk1andk1 = (k2*mu12+k1)/k1pk2                    # Cosine angle between \vec{k2}+\vec{k1} and \vec{k1}
    # Interpolate power spectra
    P1   = Plin_int(k1)
    P2   = Plin_int(k2)
    P1m2 = Plin_int(k1mk2)
    P1p2 = Plin_int(k1pk2)
    # Define terms
    term1 = 12.*F3s(k1,k1,k2,-1.+epsilon_treelevel,mu12,-mu12)             *  P1**2. * P2
    term2 = 12.*F3s(k2,k2,k1,-1.+epsilon_treelevel,mu12,-mu12)             *  P2**2. * P1
    term3 = 4.*F2s(k1mk2,k2, ang_k1mk2andk2)**2. * P2**2. * P1m2
    term4 = 4.*F2s(k1pk2,k2,-ang_k1pk2andk2)**2. * P2**2. * P1p2
    term5 = 4.*F2s(k1mk2,k1, ang_k2mk1andk1)**2. * P1**2. * P1m2
    term6 = 4.*F2s(k1pk2,k1,-ang_k2pk1andk1)**2. * P1**2. * P1p2
    term7 = 8.*F2s(k1mk2,k2, ang_k1mk2andk2)*F2s(k1mk2,k1, ang_k2mk1andk1)  *  P1 * P2 * P1m2
    term8 = 8.*F2s(k1pk2,k2,-ang_k1pk2andk2)*F2s(k1pk2,k1,-ang_k2pk1andk1)  *  P1 * P2 * P1p2
    return term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8

# ================================================================= #
# Compute multipoles of the standard tree level covariance and write them to file
# ================================================================= #

mu_range = linspace(-1.+epsilon_treelevel, 1.-epsilon_treelevel, int(Nmu_treelevel)) # for the angle average integral

cov_stdtree_mono_towrite = zeros([len(k_common), len(k_common)])
cov_stdtree_quad_towrite = zeros([len(k_common), len(k_common)])
cov_stdtree_octu_towrite = zeros([len(k_common), len(k_common)])

for i in range(len(k_common)):
    if(mod(i, int(len(k_common)/10))==0):
        print ('Up to line', i, 'out of', len(k_common),' ::  Each part gets faster though ... ')
    k2now = k_common[i]
    for j in range(i, len(k_common)):
        k1now          = k_common[j]

        integrand_mono = cov_stdtree_angle(k1now, k2now, mu_range)
        integrand_quad = cov_stdtree_angle(k1now, k2now, mu_range) * Plegendre_2(mu_range)
        integrand_octu = cov_stdtree_angle(k1now, k2now, mu_range) * Plegendre_4(mu_range)

        cov_stdtree_mono_towrite[i,j] = (1./2)*integrate.trapz(integrand_mono, mu_range) 
        cov_stdtree_quad_towrite[i,j] = (5./2)*integrate.trapz(integrand_quad, mu_range) 
        cov_stdtree_octu_towrite[i,j] = (9./2)*integrate.trapz(integrand_octu, mu_range) 

# Symmetrize for k1 <--> k2
cov_stdtree_mono_towrite = symmetrize_matrix(cov_stdtree_mono_towrite) / Volume
cov_stdtree_quad_towrite = symmetrize_matrix(cov_stdtree_quad_towrite) / Volume
cov_stdtree_octu_towrite = symmetrize_matrix(cov_stdtree_octu_towrite) / Volume

# Write
covwriter_matrix(cov_stdtree_mono_towrite, '../data_output/data_cov_stdtree_mono.dat')
covwriter_matrix(cov_stdtree_quad_towrite, '../data_output/data_cov_stdtree_quad.dat')
covwriter_matrix(cov_stdtree_octu_towrite, '../data_output/data_cov_stdtree_octu.dat')






