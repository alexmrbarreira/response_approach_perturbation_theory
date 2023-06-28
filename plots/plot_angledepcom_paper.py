import sys; sys.path.append('../'); 
from parameters import *
from functions import *

from load_data_cov import *

# ================================================================= #
# Plot map and covariance matrix slices
# ================================================================= #
labelsize = 24
ticksize  = 20
textsize  = 20
titlesize = 22
scaleaxis = 'linear'
v_min = 0.0
v_max = 0.6

fig1 = plt.figure(1, figsize=(17.5, 5.))
fig1.subplots_adjust(left=0.07, right=0.99, top=0.92, bottom=0.16, wspace = 0.25, hspace = 0.15)

fig1.add_subplot(1,3,1)
plot(k_common,  1.0e3*diagonal(cov_total_mono)/Pnl_int(k_common)**2., linestyle = 'solid'  , linewidth = 2.0, c = 'r', label = r"$\ell = 0$")
plot(k_common,  1.0e3*diagonal(cov_total_quad)/Pnl_int(k_common)**2., linestyle = 'solid'  , linewidth = 2.0, c = 'm', label = r"$\ell = 2$")
plot(k_common,  1.0e3*diagonal(cov_total_octu)/Pnl_int(k_common)**2., linestyle = 'solid'  , linewidth = 2.0, c = 'g', label = r"$\ell = 4$")
plot(k_common,  1.0e3*diagonal(cov_gauss_mono)/Pnl_int(k_common)**2., linestyle = 'dashed' , linewidth = 2.0, c = 'k', label = r"$Gaussian$")
xlabel(r'$k\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r"${\rm Cov}^{\ell}(k,k)/[P_{m}(k)]^2 \times\ 10^{3}$" , fontsize = labelsize)
xlim(0.0, 1.2)
ylim(1.0e-2, 1.0e2)
yscale('log')
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': labelsize-2}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)

fig1.add_subplot(1,3,2)
i = 15
k2_now = k_common[i]
plot(k_common,  1.0e3*(cov_total_mono[:,i]        / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"$\ell = 0$")
plot(k_common,  1.0e3*(cov_total_quad[:,i]        / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'm', label = r"$\ell = 2$")
plot(k_common,  1.0e3*(cov_total_octu[:,i]        / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'g', label = r"$\ell = 4$")
title(r"$k_2 = " + str("%.3f" % k2_now) + r"\ [h/{\rm Mpc}]$", color = 'k', fontsize = titlesize)
xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r"${\rm Cov}^{\ell}/(P_{m,1} P_{m,2}) \times\ 10^{3}$" , fontsize = labelsize)
xlim(min(k_common), max(k_common))
ylim(0.0, 1.0)
xscale('log')
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': labelsize}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)

fig1.add_subplot(1,3,3)
i = 100
k2_now = k_common[i]
plot(k_common,  1.0e3*(cov_total_mono[:,i]        / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"$\ell = 0$")
plot(k_common,  1.0e3*(cov_total_quad[:,i]        / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'm', label = r"$\ell = 2$")
plot(k_common,  1.0e3*(cov_total_octu[:,i]        / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'g', label = r"$\ell = 4$")
title(r"$k_2 = " + str("%.3f" % k2_now) + r"\ [h/{\rm Mpc}]$", color = 'k', fontsize = titlesize)
xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r"${\rm Cov}^{\ell}/(P_{m,1} P_{m,2}) \times\ 10^{3}$" , fontsize = labelsize)
xlim(min(k_common), max(k_common))
ylim(0.0, 1.0)
xscale('log')
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': labelsize}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)


fig1.savefig('fig_angledepcomp_paper.png')

show()
