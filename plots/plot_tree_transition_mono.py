import sys; sys.path.append('../'); from parameters import *
import sys; sys.path.append('../'); from functions import *

from load_data_cov import *

# ================================================================= #
# Plot
# ================================================================= #
labelsize = 20
ticksize  = 17
textsize  = 20
titlesize = 18
legendsize = 18
v_max = 0
v_min = -5

fig0 = plt.figure(0, figsize=(14., 6.))
fig0.subplots_adjust(left=0.06, right=0.99, top=0.96, bottom=0.12, wspace = 0.16, hspace = 0.32)

panel = fig0.add_subplot(1,2,1)
panel.set_aspect('equal')
title(r'${\rm log}_{10}\left[{\rm Cov}^{\rm NG, \ell = 0}_{\rm st-tree}/(P_{m,1} P_{m,2}) \times\ 10^{3}\right]$', fontsize = titlesize)
pcolor(k_common, k_common, log10(1.0e3*cov_stitchedtree_mono/(matrix_PiPj)), vmin = v_min, vmax  = v_max, cmap='jet')
colorbar().ax.tick_params(labelsize=ticksize)
plot(k_common, squeezed_transition*k_common, c = 'k', linewidth = 2, linestyle = 'dashed')
plot(k_common, (1./squeezed_transition)*k_common, c = 'k', linewidth = 2, linestyle = 'dashed')
annotate(r"$Standard\ tree$", xy = (0.45, 0.60), xycoords='axes fraction', rotation = 00, xytext = (0.45, 0.60), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize+2)
annotate(r"$Response\ tree$", xy = (0.03, 0.90), xycoords='axes fraction', rotation = 00, xytext = (0.03, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize+2)
annotate(r"$Response\ tree$", xy = (0.60, 0.12), xycoords='axes fraction', rotation = 00, xytext = (0.60, 0.12), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize+2)
xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$k_2\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
xlim(min((k_common)), max((k_common)))
ylim(min((k_common)), max((k_common)))
#xscale('log')
#yscale('log')
xticks(size = ticksize)
yticks(size = ticksize)

fig0.add_subplot(2,2,2)
i = 15
axvline(squeezed_transition*k_common[i], c = 'k', linewidth = 2, linestyle = 'dashed')
axvline((1./squeezed_transition)*k_common[i], c = 'k', linewidth = 2, linestyle = 'dashed')
plot(k_common,  1.0e3*(cov_stdtree_mono[:,i]      / Pnl_int(k_common)/ Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'b', label = r"$Standard\ tree$")
plot(k_common,  1.0e3*(cov_sqtree_mono[:,i]       / Pnl_int(k_common)/ Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'g', label = r"$Response\ tree$")
plot(k_common,  1.0e3*(cov_stitchedtree_mono[:,i] / Pnl_int(k_common)/ Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"$Stitched\ tree$")
#plot(k_common,  1.0e3*(cov_stdtree_mono[:,i]      /Plin_int(k_common)/Plin_int(k_common[i])), linestyle = 'dashed', linewidth = 1.5, c = 'm', label = r"$Standard\ with\ corr.$")

#annotate(r"$k_2 = "+str(k_common[i])+r"\ h/{\rm Mpc}$", xy = (0.04, 0.85), xycoords='axes fraction', color = 'k', fontsize = textsize)
annotate(r"$k_2 = 0.086"+r"\ h/{\rm Mpc}$", xy = (0.02, 0.85), xycoords='axes fraction', color = 'k', fontsize = textsize)
xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
#ylabel(r'${\rm Cov}^{\rm NG,\ell=0}_{12}/P_{m,1}/P_{m,2}\ \times\ 10^{3}$'        , fontsize = labelsize)
xlim(min((k_common)), max((k_common)))
ylim(0.0, 1.0)
xscale('log')
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legendsize}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)

fig0.add_subplot(2,2,4)
i = 100
axvline(squeezed_transition*k_common[i], c = 'k', linewidth = 2, linestyle = 'dashed')
axvline((1./squeezed_transition)*k_common[i], c = 'k', linewidth = 2, linestyle = 'dashed')
plot(k_common,  1.0e3*(cov_stdtree_mono[:,i]      / Pnl_int(k_common)/ Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'b', label = r"$Standard\ tree$")
plot(k_common,  1.0e3*(cov_sqtree_mono[:,i]       / Pnl_int(k_common)/ Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'g', label = r"$Response\ tree$")
plot(k_common,  1.0e3*(cov_stitchedtree_mono[:,i] / Pnl_int(k_common)/ Pnl_int(k_common[i])), linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"$Stitched\ tree$")
#plot(k_common,  1.0e3*(cov_stdtree_mono[:,i]      /Plin_int(k_common)/Plin_int(k_common[i])), linestyle = 'dashed', linewidth = 1.5, c = 'm', label = r"$Standard\ with\ corr.$")

#annotate(r"$k_2 = "+str(k_common[i])+r"\ h/{\rm Mpc}$", xy = (0.04, 0.85), xycoords='axes fraction', color = 'k', fontsize = textsize)
annotate(r"$k_2 = 0.49"+r"\ h/{\rm Mpc}$", xy = (0.02, 0.85), xycoords='axes fraction', color = 'k', fontsize = textsize)
xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ {\rm Cov}^{\rm NG, \ell=0}_{\rm st-tree}/(P_{m,1} P_{m,2})\ \times\ 10^{3}$'        , fontsize = labelsize)
xlim(min((k_common)), max((k_common)))
ylim(0.0, 1.0)
xscale('log')
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legendsize}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)

fig0.savefig('fig_tree_transition_mono.png')

show()
