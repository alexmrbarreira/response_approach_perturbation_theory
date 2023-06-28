import sys; sys.path.append('../'); 
from parameters import *
from functions import *

from load_data_cov import *

# ================================================================= #
# Plot Cov(k1, k3) for a few k3 values monopole
# ================================================================= #
fig0 = plt.figure(0, figsize=(17.5, 10.))
fig0.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.10, wspace = 0.19, hspace = 0.38)

i_list = [7,8,9,10,11,12] #Chose indices of k_common here by hand
counter=0
for i in i_list: 
    fig0.add_subplot(2,3,counter+1)

    plot(k_common,  1.0e3*(cov_sqtree_mono[:,i] /Pnl_int(k_common)/Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"${\rm Cov}^{\rm NG, sq}_{\ell=0}(k_1, k_2)$")

    x_fill = linspace((1./squeezed_transition)*k_common[i], squeezed_transition*k_common[i], 5)
    fill_between(x_fill, -2.0e0, 2.0e0, color = 'grey', alpha = 0.7)
    plot(x_fill, 1.0e99*ones(len(x_fill)), color = 'grey', linewidth = 10., label = r'$\approx\ non-squeezed$')

    title(r"$k_2 = " + str("%.3f" % k_common[i]) + r"\ [h/{\rm Mpc}]$", color = 'k', fontsize = 18)
    xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = 18)
    ylabel(r"${\rm Cov}_{\ell=0}(k_1, k_2)/P_m(k_1)/P_m(k_2)\ \times\ 10^3$" , fontsize = 18)
    xlim(min(k_common), max(k_common))
    ylim(0.0, 1.0)
    xscale('log')
    xticks(size = 18)
    yticks(size = 18)
    params = {'legend.fontsize': 16}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)
    counter += 1

# ================================================================= #
# Plot Cov(k1, k3) for a few k3 values quadrupole
# ================================================================= #
fig1 = plt.figure(1, figsize=(17.5, 5.))
fig1.subplots_adjust(left=0.06, right=0.99, top=0.92, bottom=0.17, wspace = 0.20, hspace = 0.38)
i_list = [7,9,12] #Chose indices of k_common here by hand
i_list = [7,8,9] #Chose indices of k_common here by hand
counter = 0
for i in i_list:
    fig1.add_subplot(1,3,counter+1)

    plot(k_common,  1.0e3*(cov_sqtree_mono[:,i] /Pnl_int(k_common)/Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"$\ell=0$")
    plot(k_common,  1.0e3*(cov_sqtree_quad[:,i] /Pnl_int(k_common)/Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'g', label = r"$\ell=2$")
    plot(k_common,  1.0e3*(cov_sqtree_octu[:,i] /Pnl_int(k_common)/Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'c', label = r"$\ell=4$")

    x_fill = linspace((1./squeezed_transition)*k_common[i], squeezed_transition*k_common[i], 5)
    fill_between(x_fill, -2.0e0, 2.0e0, color = 'grey', alpha = 0.7)
    plot(x_fill, 1.0e99*ones(len(x_fill)), color = 'grey', linewidth = 10., label = r'$\approx\ non-squeezed$')

    title(r"$k_2 = " + str("%.3f" % k_common[i]) + r"\ [h/{\rm Mpc}]$", color = 'k', fontsize = 18)
    xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = 18)
    ylabel(r"${\rm Cov}^{\rm NG, sq}_{\ell}(k_1, k_2)/P_m(k_1)/P_m(k_2)\ \times\ 10^{3}$" , fontsize = 18)
    xlim(min(k_common), max(k_common))
    ylim(0.0, 1.0e0)
    xscale('log')
    xticks(size = 18)
    yticks(size = 18)
    params = {'legend.fontsize': 17}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 2)
    counter += 1


fig0.savefig('fig_sqtreecov_mono_someslices.png')

fig1.savefig('fig_sqtreecov_multi_someslices.png')

#show()
