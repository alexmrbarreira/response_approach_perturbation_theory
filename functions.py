from parameters import *

# ==================================================================== #
# Some useful functions 
# ==================================================================== #
def Spline3(xx,yy,Order,X):
    ff=interpolate.splrep(xx,yy,k=3)
    ss=interpolate.splev(X,ff,Order)
    return ss

def finidiff_2nd(xx,yy): ## Assumes uniform spacing of xx || Computes the second derivative with finite differences
    h = xx[1] - xx[0]
    d2yy = zeros(len(xx))
    # Do Spline method for first two entries
    d2yy[0] = Spline3(xx,yy,2,xx[0])
    d2yy[1] = Spline3(xx,yy,2,xx[1])
    # Do central finite differences
    index_now = array(range(2, len(xx)-2))
    d2yy[index_now] = ( (-1./12)*yy[index_now-2] + (4./3)*yy[index_now-1] + (-5./2)*yy[index_now] + (4./3)*yy[index_now+1] + (-1./12)*yy[index_now+2] ) / h**2.
    # Do Spline method for last two entries
    d2yy[-2] = Spline3(xx,yy,2,xx[-2])
    d2yy[-1] = Spline3(xx,yy,2,xx[-1])
    return d2yy

# Function to symmetrize calculation of the covariance matrices: speed up things by computing only the triangle k1>=k2
# Useful for squeezed covariance because the matrix has to be symmetric wrt k1=k2, but things depends on whether k1><k2; we compute the matrix assuming k1>k2, but then symmetrize it 
def symmetrize_matrix(matrix):
    if(len(matrix[:,0]) != len(matrix[0,:])):
        print ('Not square matrix!'); quit()
    n = len(matrix[:,0])
    new_matrix = copy(matrix)
    for j in range(n):
        for i in range(1+j, n):
            new_matrix[i, j] = matrix[j, i]
    return new_matrix

# ==================================================================== #
# Create interpolators for linear power spectra
# ==================================================================== #
klin      = f_pk_camb[:,0]
Plin      = f_pk_camb[:,1]
Plin_int  = interpolate.interp1d(klin,   Plin)
# Make 1st derivative interpolator (Spline function suffices)
dPlin  = zeros(len(klin)) #This will be dP/dk
for i in range(len(klin)):
    dPlin[i]  = Spline3(klin,Plin,1,klin[i])
dPlin_int  = interpolate.interp1d(klin,  dPlin)

# Make 2nd derivative interpolator (Spline function suffices)
ddPlin  = zeros(len(klin)) #This will be d^2P/dk^2
for i in range(len(klin)):
    ddPlin[i]  = Spline3(klin,dPlin,1,klin[i])

# Smooth out numerical noise
def smoother(ddP, kk, sigma, kmin, kmax):
    smoothed =  scipy.ndimage.filters.gaussian_filter(ddP, sigma)
    out = copy(ddP)
    out[where(kk<kmin)] = smoothed[where(kk<kmin)]
    out[where(kk>kmax)] = smoothed[where(kk>kmax)]
    return out
ddPlin_smt = smoother(ddPlin, klin, 1.5, 0.05, 0.6)
ddPlin_int  = interpolate.interp1d(klin,  ddPlin_smt)

# ==================================================================== #
# Load nonlinear P(k) file from Coyote and create interpolator
# ==================================================================== #

knl      = f_pk_coyo[:,0]/h_hubble
Pnl      = f_pk_coyo[:,1]*h_hubble**3.
Pnl_int  = interpolate.interp1d(knl,   Pnl)

# Make 1st derivative interpolator (Spline function suffices)
dPnl  = zeros(len(knl)) #This will be dP/dk
for i in range(len(knl)):
    dPnl[i]  = Spline3(knl,Pnl,1,knl[i])
dPnl_int  = interpolate.interp1d(knl,  dPnl)

# Make 2nd derivative interpolator (use explicit finite difference function)
knl_fini  = exp(linspace(log(min(knl)+tinycorr), log(max(knl)), len(knl)))
Pnl_fini  = Pnl_int(knl_fini)
ddPnl     = finidiff_2nd(log(knl_fini), Pnl_fini)            # This will be P'', with ' = d/dk
ddPnl     = (ddPnl/knl_fini - dPnl_int(knl_fini))/knl_fini
ddPnl_int = interpolate.interp1d(knl_fini, ddPnl)

# ==================================================================== #
# Define tree-level response functions  (Lagrangian and Eulerian -- indicated in the functions names as l and e)
# ==================================================================== #

def fk1(k): # This is kP'/P
    return k     *  dPlin_int(k) / Plin_int(k)
def fk2(k): # This is k^2 P''/P
    return k**2. * ddPlin_int(k) / Plin_int(k)
def R1l_tree(k): #This is R_1(k)
    return 47./21 - (1./3)*fk1(k)
def R2l_tree(k): #This is R_K(k)
    return 8./7 - fk1(k)
def R3l_tree(k): #This is R_2(k)
    return 8420./1323 - (100./63)*fk1(k) + (1./9)*fk2(k)
def R4l_tree(k): #This is R_Kdelta(k)
    return 1348./441 - (55./21)*fk1(k) + (1./3)*fk2(k)
def R5l_tree(k): #This is R_K^2(k)
    return 20./63 + (1./14)*fk1(k)
def R6l_tree(k): #This is R_K.K(k)
    return -20./21 + (1./2)*fk1(k)
def R7l_tree(k): #This is R_KK(k)
    return 328./147 - (23./14)*fk1(k) + (1./2)*fk2(k)
def R8l_tree(k): #This is R_td(k)
    return 8./63 - (1./7)*fk1(k)
def R1e_tree(k): #This is R_1(k)
    return (47./21) - (1./3)*fk1(k) 
def R2e_tree(k): #This is R_K(k)
    return (8./7) - fk1(k)
def R3e_tree(k): #This is R_2(k)
    return (74./27) - (22./21)*fk1(k) + (1./9)*fk2(k)
def R4e_tree(k): #This is R_Kdelta(k)
    return (1012./441) - (41./21)*fk1(k) + (1./3)*fk2(k)
def R5e_tree(k): #This is R_K^2(k)
    return (26./441) - (1./6)*fk1(k)
def R6e_tree(k): #This is R_K.K(k)
    return (-44./21) + (3./2)*fk1(k)
def R7e_tree(k): #This is R_KK(k)
    return (328./147) - (23./14)*fk1(k) + (1./2)*fk2(k)
def R8e_tree(k): #This is R_td(k)
    return (-184./441) + (1./3)*fk1(k)

# ==================================================================== #
# Define extrapolated response functions 
# ==================================================================== #
def fk1nl(k): # This is kP'/P
    return k     *  dPnl_int(k) / Pnl_int(k)
def fk2nl(k): # This is k^2 P''/P
    return k**2. * ddPnl_int(k) / Pnl_int(k)

# 1) Define Lagrangian ones from simulations and those which are the same as Eulerian
def R1l_nonl(k): #This is R_1(k)
    return R1_rnl_int(k)
def R2l_nonl(k): #This is R_K(k)
    return ((12./13)*G1_gnl_int(k) - fk1nl(k))
def R3l_nonl(k): #This is R_2(k)
    return R2_rnl_int(k)
def R7l_nonl(k): #This is R_KK(k)
    return (738./2105)*(34./21 + 2.*G1_gnl_int(k) + G2_gnl_int(k)) - (207./200)*(16./21 + (2./3)*G1_gnl_int(k))*fk1nl(k) + (1./2)*fk2nl(k)

# 2) Define Eulerian ones using the extrapolation as described in paper
def R1e_nonl(k): #This is R_1(k)
    return R1l_nonl(k)
def R2e_nonl(k): #This is R_K(k)
    return ((12./13)*G1_gnl_int(k) - fk1nl(k))
def R3e_nonl(k): #This is R_2(k)
    return R3l_nonl(k) - (34./21)*R1e_nonl(k)
def R4e_nonl(k): #This is R_Kdelta(k)
    return (1518./1813)*((8./21)*G1_gnl_int(k) + G2_gnl_int(k)) + (41./22)*(-2./9 - (2./3)*G1_gnl_int(k))*fk1nl(k) + (1./3)*fk2nl(k)
def R5e_nonl(k): #This is R_K^2(k)
    return (1./21)*G1_gnl_int(k) - (1./6)*fk1nl(k)
def R6e_nonl(k): #This is R_K.K(k)
    return -(22./13)*G1_gnl_int(k) + (3./2)*fk1nl(k)
def R7e_nonl(k): #This is R_KK(k)
    return (1476./1813)*((8./21)*G1_gnl_int(k) + G2_gnl_int(k)) + (69./44)*(-2./9 - (2./3)*G1_gnl_int(k))*fk1nl(k) + (1./2)*fk2nl(k)
def R8e_nonl(k): #This is R_td(k)
    return -(92./273)*G1_gnl_int(k) + (1./3)*fk1nl(k)

# 3) Define remaining Lagrangian ones from the relations to the Eulerian ones
def R4l_nonl(k): #This is R_Kdelta(k)
    return (2./3)*R2e_nonl(k) + R4e_nonl(k)
def R5l_nonl(k): #This is R_K^2(k)
    return (2./7)*R1e_nonl(k) + (-1./3)*R2e_nonl(k) + R5e_nonl(k)
def R6l_nonl(k): #This is R_K.K(k)
    return R2e_nonl(k) + R6e_nonl(k)
def R8l_nonl(k): #This is R_td(k)
    return (10./21)*R2e_nonl(k) + R8e_nonl(k)

# ==================================================================== #
# Define total second order response functions (Lagrangian and Eulerian)
# ==================================================================== #

def curlyR2l_tree(k, m1, m2, m12): #\mathcal{R}_2 at tree level
    term1 = m12
    term2 = (1./2)*m12*(m1**2 + m2**2. - 2./3)
    term3 = 1./2
    term4 = (1./2)*(m1**2 + m2**2. - 2./3)
    term5 = m12**2. - 1./3
    term6 = m1*m2*m12 - (1./3)*m1**2. - (1./3)*m2**2. + 1./9
    term7 = m1**2.*m2**2. - (1./3)*(m1**2. + m2**2.) + 1./9
    term8 = (3./2)*((1./2)*(1-m12)*(m1 + m2)**2. - (1./3)*(1-m12**2.))
    return R1l_tree(k)*term1 + R2l_tree(k)*term2 + R3l_tree(k)*term3 + R4l_tree(k)*term4 + R5l_tree(k)*term5 + R6l_tree(k)*term6 + R7l_tree(k)*term7 + R8l_tree(k)*term8 

def curlyR2l_nonl(k, m1, m2, m12): #\mathcal{R}_2 with nonlinear extrapolation of coefficients
    term1 = m12
    term2 = (1./2)*m12*(m1**2 + m2**2. - 2./3)
    term3 = 1./2
    term4 = (1./2)*(m1**2 + m2**2. - 2./3)
    term5 = m12**2. - 1./3
    term6 = m1*m2*m12 - (1./3)*m1**2. - (1./3)*m2**2. + 1./9
    term7 = m1**2.*m2**2. - (1./3)*(m1**2. + m2**2.) + 1./9
    term8 = (3./2)*((1./2)*(1-m12)*(m1 + m2)**2. - (1./3)*(1-m12**2.))
    return R1l_nonl(k)*term1 + R2l_nonl(k)*term2 + R3l_nonl(k)*term3 + R4l_nonl(k)*term4 + R5l_nonl(k)*term5 + R6l_nonl(k)*term6 + R7l_nonl(k)*term7 + R8l_nonl(k)*term8

def curlyR2e_tree(k, m1, m2, m12): #\mathcal{R}_2 at tree level
    term1 = 5./7 + m12 + (2./7)*m12**2.
    term2 = m1*m2*m12 - (1./3)*m12**2. + (5./14)*(1-m12)*(m1+m2)**2. - (5./21)*(1-m12**2.) + (1./2)*m12*(m1**2 + m2**2. - 2./3)
    term3 = 1./2
    term4 = (1./2)*(m1**2 + m2**2. - 2./3)
    term5 = m12**2. - 1./3
    term6 = m1*m2*m12 - (1./3)*m1**2. - (1./3)*m2**2. + 1./9
    term7 = m1**2.*m2**2. - (1./3)*(m1**2. + m2**2.) + 1./9
    term8 = (3./2)*((1./2)*(1-m12)*(m1 + m2)**2. - (1./3)*(1-m12**2.))
    return R1e_tree(k)*term1 + R2e_tree(k)*term2 + R3e_tree(k)*term3 + R4e_tree(k)*term4 + R5e_tree(k)*term5 + R6e_tree(k)*term6 + R7e_tree(k)*term7 + R8e_tree(k)*term8

def curlyR2e_nonl(k, m1, m2, m12): #\mathcal{R}_2 with nonlinear extrapolation of coefficients
    term1 = 5./7 + m12 + (2./7)*m12**2.
    term2 = m1*m2*m12 - (1./3)*m12**2. + (5./14)*(1-m12)*(m1+m2)**2. - (5./21)*(1-m12**2.) + (1./2)*m12*(m1**2 + m2**2. - 2./3)
    term3 = 1./2
    term4 = (1./2)*(m1**2 + m2**2. - 2./3)
    term5 = m12**2. - 1./3
    term6 = m1*m2*m12 - (1./3)*m1**2. - (1./3)*m2**2. + 1./9
    term7 = m1**2.*m2**2. - (1./3)*(m1**2. + m2**2.) + 1./9
    term8 = (3./2)*((1./2)*(1-m12)*(m1 + m2)**2. - (1./3)*(1-m12**2.))
    return R1e_nonl(k)*term1 + R2e_nonl(k)*term2 + R3e_nonl(k)*term3 + R4e_nonl(k)*term4 + R5e_nonl(k)*term5 + R6e_nonl(k)*term6 + R7e_nonl(k)*term7 + R8e_nonl(k)*term8

# ==================================================================== #
# Define alpha and beta mode coupling functions 
# ==================================================================== #

def alpha(q,k,mu):
    return 1. + k*mu/q

def beta(q,k,mu):
    return (mu/2.)*(k/q + q/k) + mu**2.

# ==================================================================== #
# Define F2 and G2 kernels and symmetrize them
# ==================================================================== #

def F2(q1,q2,mu12):
    return (1./7)*(5.*alpha(q1,q2,mu12) + 2.*beta(q1,q2,mu12))

def G2(q1,q2,mu12):
    return (1./7)*(3.*alpha(q1,q2,mu12) + 4.*beta(q1,q2,mu12))

def F2s(q1,q2,mu12):
    return (1./2)*(F2(q1,q2,mu12) + F2(q2,q1,mu12))

def G2s(q1,q2,mu12):
    return (1./2)*(G2(q1,q2,mu12) + G2(q2,q1,mu12))

# ==================================================================== #
# Define F3 and symmetrize it
# ==================================================================== #

def F3(q1,q2,q3,mu12,mu13,mu23):
    term1 = 2.* beta(q1,sqrt(q2**2.+q3**2.+2.*q2*q3*mu23),(q2*mu12+q3*mu13)/sqrt(q2**2.+q3**2.+2.*q2*q3*mu23)) * G2(q2,q3,mu23)
    term2 = 7.*alpha(q1,sqrt(q2**2.+q3**2.+2.*q2*q3*mu23),(q2*mu12+q3*mu13)/sqrt(q2**2.+q3**2.+2.*q2*q3*mu23)) * F2(q2,q3,mu23)
    term3 = 2.* beta(sqrt(q1**2.+q2**2.+2.*q1*q2*mu12),q3,(q1*mu13+q2*mu23)/sqrt(q1**2.+q2**2.+2.*q1*q2*mu12)) * G2(q1,q2,mu12)
    term4 = 7.*alpha(sqrt(q1**2.+q2**2.+2.*q1*q2*mu12),q3,(q1*mu13+q2*mu23)/sqrt(q1**2.+q2**2.+2.*q1*q2*mu12)) * G2(q1,q2,mu12)
    return (1./18)*(term1 + term2 + term3 + term4)

def F3s(q1,q2,q3,mu12,mu13,mu23):
    out = F3(q1,q2,q3,mu12,mu13,mu23) + F3(q1,q3,q2,mu13,mu12,mu23) + F3(q2,q1,q3,mu12,mu23,mu13) + F3(q2,q3,q1,mu23,mu12,mu13) + F3(q3,q1,q2,mu13,mu23,mu12) + F3(q3,q2,q1,mu23,mu13,mu12)
    return out/6.

# ==================================================================== #
# Functions to write to files 
# ==================================================================== #

def covwriter_matrix(cov_2d_array, filepath):
    fout = open(filepath, 'w')
    leny = len(cov_2d_array[:,0])
    lenx = len(cov_2d_array[0,:])
    for i in range(leny):
        for j in range(lenx):
            fout.write(str(cov_2d_array[i,j])); fout.write(' ');
        fout.write('\n')
    fout.close()
    return 0

# ==================================================================== #
# A few Legendre polynomials (even only)
# ==================================================================== #

def Plegendre_0(x):
    return 1.
def Plegendre_2(x):
    return (3.*x**2. - 1.)/2.
def Plegendre_4(x):
    return (35.*x**4. - 30.*x**2. + 3.)/8
def Plegendre_6(x):
    return (231.*x**6. - 315.*x**4. + 105.*x**2. - 5)/16.
def Plegendre_8(x):
    return (6435.*x**8. - 12012.*x**6. + 6930.*x**4. - 1260.*x**2. + 35.)/128.
def Plegendre_10(x):
    return (46189.*x**10 - 109395.*x**8. + 90090.*x**6. - 30030.*x**4. + 3465.*x**2. - 63.)/256. 

# ==================================================================== #
# Define common k-array to use in interpolations (define this always at the end of this file to have all used k arrays)
# ==================================================================== #
k_sims       = load(localpath + 'data_input/k_sims.npy')
k_common_min = max(array([min(k_gnl), min(k_rnl), min(klin), min(knl)]))
k_common_max = min(array([max(k_gnl), max(k_rnl), max(klin), max(knl)]))
k_common     = k_sims[where((k_sims > k_common_min) & (k_sims < k_common_max))]


