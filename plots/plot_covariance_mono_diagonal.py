import sys; sys.path.append('../'); from parameters import *
import sys; sys.path.append('../'); from functions import *

from load_data_cov import *

# ================================================================= #
# Plot
# ================================================================= #
labelsize = 18
ticksize  = 15
textsize  = 18
titlesize = 16
legendsize = 15
v_max = 1.0e0
v_min = 1.0e-5

fig0 = plt.figure(0, figsize=(8., 5.))
fig0.subplots_adjust(left=0.16, right=0.97, top=0.97, bottom=0.15, hspace = 0.15)#, space = 0.2)

fig0.add_subplot(1,1,1)
plot(k_common,  1.0e3*diagonal(cov_gauss_mono        / matrix_PiPj), linestyle = 'solid' , linewidth = 2.0, c = 'c', label = r"$Gaussian$")
plot(k_common,  1.0e3*diagonal(cov_stitchedtree_mono / matrix_PiPj), linestyle = 'solid' , linewidth = 2.0, c = 'b', label = r"$Stitched\ tree$")
plot(k_common,  1.0e3*diagonal(cov_1loop_mono        / matrix_PiPj), linestyle = 'solid' , linewidth = 2.0, c = 'g', label = r"$Response\ 1-loop$")
plot(k_common,  1.0e3*diagonal(cov_total_mono        / matrix_PiPj), linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"$Total$")
annotate(r"$z=0$", xy = (0.55,0.60), xycoords='axes fraction', color = 'k', fontsize = 24)
xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'${\rm Cov}^{\ell = 0}(k,k)/[P_{m}(k)]^2 \times\ 10^{3}$'        , fontsize = labelsize)
xlim(min((k_common)), max((k_common)))
ylim(1.0e-3, 1.0e4)
xscale('log')
yscale('log')
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legendsize+4}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 2)

fig0.savefig('fig_covariance_diagonal_mono.png')

show()
