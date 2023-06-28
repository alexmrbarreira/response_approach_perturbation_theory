import sys; sys.path.append('../'); 
from parameters import *
from functions import *

from load_data_cov import *

# ================================================================= #
# Plot Covariance matrices 
# ================================================================= #
labelsize = 24
ticksize  = 22
textsize  = 28
titlesize = 28
v_min = 0.0
v_max = 0.6
scaleaxis = 'linear'

v_min = -0.3
v_max = 1.0
fig1 = plt.figure(1, figsize=(16., 16.))
fig1.subplots_adjust(left=0.08, right=0.985, top=0.98, bottom=0.04, wspace = 0.15, hspace = 0.10)

panel = fig1.add_subplot(2,2,1)
panel.set_aspect('equal')
title(r'$r_{\ell=0}(k_1, k_2)$', fontsize = titlesize)
pcolor(k_common, k_common, cov_stitchedtree_mono/sqrt(matrix_CiiCjj_mono), vmin = v_min, vmax = v_max, cmap='jet')
colorbar().ax.tick_params(labelsize=ticksize-2)
annotate(r"$Stitched\ tree$", xy = (0.04, 0.90), xycoords='axes fraction',
xytext = (0.04, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize)
xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$k_2\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
xlim(min(k_common), max(k_common))
ylim(min(k_common), max(k_common))
xscale(scaleaxis)
yscale(scaleaxis)
xticks(size = ticksize)
yticks(size = ticksize)

panel = fig1.add_subplot(2,2,2)
panel.set_aspect('equal')
title(r'$r_{\ell=0}(k_1, k_2)$', fontsize = titlesize)
pcolor(k_common, k_common, cov_1loop_mono/sqrt(matrix_CiiCjj_mono), vmin = v_min, vmax = v_max, cmap='jet')
colorbar().ax.tick_params(labelsize=ticksize-2)
annotate(r"$Response\ 1-loop$", xy = (0.04, 0.90), xycoords='axes fraction',
xytext = (0.04, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize)
xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$k_2\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
xlim(min(k_common), max(k_common))
ylim(min(k_common), max(k_common))
xscale(scaleaxis)
yscale(scaleaxis)
xticks(size = ticksize)
yticks(size = ticksize)

panel = fig1.add_subplot(2,2,3)
panel.set_aspect('equal')
title(r'$r_{\ell=0}(k_1, k_2)$', fontsize = titlesize)
pcolor(k_common, k_common, cov_total_mono/sqrt(matrix_CiiCjj_mono), vmin = v_min, vmax = v_max, cmap='jet')
colorbar().ax.tick_params(labelsize=ticksize-2)
annotate(r"$Total$", xy = (0.04, 0.90), xycoords='axes fraction',
xytext = (0.04, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize)
xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$k_2\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
xlim(min(k_common), max(k_common))
ylim(min(k_common), max(k_common))
xscale(scaleaxis)
yscale(scaleaxis)
xticks(size = ticksize)
yticks(size = ticksize)

fig1.savefig('fig_covariance_maps_corr_mono.png')

show()
