import sys; sys.path.append('../'); 
from parameters import *
from functions import *

from load_data_cov import *

# ================================================================= #
# Plot all slices 
# ================================================================= #
labelsize = 20
ticksize  = 20
textsize  = 20
titlesize = 20

indstouse = [6, 9, 13, 19, 28, 39, 62, 88, 206] #Select, by hand, indices of k_common array for the respective slices

fig0 = plt.figure(0, figsize=(17.5, 20.))
fig0.subplots_adjust(left=0.07, right=0.99, top=0.98, bottom=0.04, wspace = 0.23, hspace = 0.25)
subplotindex = 1
for i in indstouse:

    k2_now = k_common[i]

    fig0.add_subplot(3,3,subplotindex)
    #plot(k_common,  1.0e3*(cov_gauss_mono[:,i]        / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'c', label = r"$Gaussian$")
    plot(k_common,  1.0e3*(cov_stitchedtree_mono[:,i] / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'b', label = r"$Stitched\ tree$")
    plot(k_common,  1.0e3*(cov_1loop_mono[:,i]        / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'g', label = r"$Response\ 1-loop$")
    plot(k_common,  1.0e3*(cov_total_mono[:,i]        / Pnl_int(k_common)/Pnl_int(k2_now)), linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"$Total$")

    title(r"$k_2 = " + str("%.3f" % k2_now) + r"\ [h/{\rm Mpc}]$", color = 'k', fontsize = titlesize)
    xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
    ylabel(r"${\rm Cov}^{\ell = 0}/(P_{m,1} P_{m,2}) \times\ 10^{3}$" , fontsize = labelsize)
    xlim(min(k_common), max(k_common))
    ylim(0.0, 1.0)
    xscale('log')
    xticks(size = ticksize)
    yticks(size = ticksize)
    params = {'legend.fontsize': labelsize}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)

    subplotindex += 1

fig0.savefig('fig_covariance_slices_mono.png')

show()
