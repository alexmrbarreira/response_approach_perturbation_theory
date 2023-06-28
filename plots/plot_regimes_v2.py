import sys; sys.path.append('../'); 
from parameters import *
from functions import *

from load_data_cov import *

squeezeness_1 = 0.25
squeezeness_2 = 0.65
ktoplot_1     = k_common[where(k_common > min(k_common)/squeezeness_1)]
ktoplot_2     = k_common[where(k_common > min(k_common)/squeezeness_2)]

cov_total_mono -= cov_gauss_mono

cov_total_mono_int        = interpolate.RectBivariateSpline(k_common, k_common, cov_total_mono       ,  kx=2, ky=2, s=0)
cov_1loop_mono_int        = interpolate.RectBivariateSpline(k_common, k_common, cov_1loop_mono       ,  kx=2, ky=2, s=0)
cov_stitchedtree_mono_int = interpolate.RectBivariateSpline(k_common, k_common, cov_stitchedtree_mono,  kx=2, ky=2, s=0)

# ================================================================= #
# Plot
# ================================================================= #
labelsize = 20
ticksize  = 17
textsize  = 19
titlesize = 18
legendsize = 15
#v_max = 0.3
#v_min = 0.05
v_max = 0
v_min = -1.3

fig0 = plt.figure(0, figsize=(14., 6.))
fig0.subplots_adjust(left=0.06, right=0.985, top=0.97, bottom=0.12, wspace = 0.16, hspace = 0.32)

panel = fig0.add_subplot(1,2,1)
panel.set_aspect('equal')
title(r'${\rm log}_{10}\left[{\rm Cov}^{\rm NG, \ell = 0}/(P_{m,1} P_{m,2}) \times\ 10^{3}\right]$', fontsize = titlesize)
pcolor(k_common, k_common, log10(1.0e3*cov_total_mono/(matrix_PiPj)), vmin = v_min, vmax  = v_max, cmap='jet')
colorbar().ax.tick_params(labelsize=ticksize)

plot(ktoplot_1, squeezeness_1*ktoplot_1, linestyle = 'dotted', linewidth = 2, c = 'k')
plot(ktoplot_2, squeezeness_2*ktoplot_2, linestyle = 'dotted', linewidth = 2, c = 'k')

kfine = 10.**linspace(-4, 1, 10000)
## Draw line of regime II
plot(kfine[where(kfine<=squeezed_transition*0.075)], (1./squeezed_transition)*kfine[where(kfine<=squeezed_transition*0.075)] , linestyle = 'solid', c = 'k', linewidth = 2.0)
plot(kfine[where(kfine>=squeezed_transition*0.075)], 0.075*ones(len(kfine[where(kfine>=squeezed_transition*0.075)]))         , linestyle = 'solid', c = 'k', linewidth = 2.0)
plot(kfine[where(kfine<=0.075)], 2.0*kfine[where(kfine<=0.075)]          , linestyle = 'solid', c = 'k', linewidth = 2.0)
plot([0.075, 0.075], [2.0*0.075, 10.]                                     , linestyle = 'solid', c = 'k', linewidth = 2.0)

## Draw line of regime III
plot(kfine[where(kfine<=0.3*squeezed_transition)], (1./squeezed_transition)*kfine[where(kfine<=0.3*squeezed_transition)]  , linestyle = 'solid', c = 'k', linewidth = 2.0)
plot(kfine[where(kfine<=0.3)], squeezed_transition    *kfine[where(kfine<=0.30)] , linestyle = 'solid', c = 'k', linewidth = 2.0)
plot([0.30*squeezed_transition, kfine[-1]], [0.30, 0.30], linestyle = 'solid', c = 'k', linewidth = 2.0)
plot([0.30, 0.30], [0.30*squeezed_transition, kfine[-1]], linestyle = 'solid', c = 'k', linewidth = 2.0)

## Draw line of regime I
plot([0.075, squeezed_transition*0.075], [0.075, 0.075],linestyle = 'solid', c = 'k', linewidth = 2.0)
plot([0.075, 0.075], [0.075, squeezed_transition*0.075],linestyle = 'solid', c = 'k', linewidth = 2.0)

## Draw line of regime IV
plot([0.30, squeezed_transition*0.30], [0.30, 0.30],linestyle = 'solid', c = 'k', linewidth = 2.0)
plot([0.30, 0.30], [0.30, squeezed_transition*0.30],linestyle = 'solid', c = 'k', linewidth = 2.0)

annotate(r"$I$",xy=(0.06,0.06),xycoords='axes fraction',rotation=0,xytext=(0.06,0.06),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)

annotate(r"$II$",xy=(0.50,0.035),xycoords='axes fraction',rotation=0,xytext=(0.50,0.035),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)
annotate(r"$II$",xy=(0.035,0.50),xycoords='axes fraction',rotation=0,xytext=(0.035,0.50),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)

annotate(r"$III$",xy=(0.75,0.25),xycoords='axes fraction',rotation=0,xytext=(0.75,0.25),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)
annotate(r"$III$",xy=(0.25,0.75),xycoords='axes fraction',rotation=0,xytext=(0.25,0.75),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)

annotate(r"$V$",xy=(0.65,0.65),xycoords='axes fraction',rotation=0,xytext=(0.65,0.65),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)

annotate(r"$IV$",xy=(0.30,0.30),xycoords='axes fraction',rotation=0,xytext=(0.30,0.30),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)

xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$k_2\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
#xlim(0.*min((k_common)), max((k_common)))
#ylim(0.*min((k_common)), max((k_common)))
xlim(0., 0.7)
ylim(0., 0.7)
#xscale('log')
#yscale('log')
xticks(size = ticksize)
yticks(size = ticksize)

fig0.add_subplot(2,2,4)

stitchedtree_curve = zeros(len(ktoplot_1))
oloop_curve        = zeros(len(ktoplot_1))
total_curve        = zeros(len(ktoplot_1))
for i in range(len(ktoplot_1)):
    stitchedtree_curve[i] = 1.0e3*(cov_stitchedtree_mono_int(ktoplot_1[i],squeezeness_1*ktoplot_1[i]) / Pnl_int(ktoplot_1[i]) / Pnl_int(squeezeness_1*ktoplot_1[i]))
    oloop_curve[i]        = 1.0e3*(       cov_1loop_mono_int(ktoplot_1[i],squeezeness_1*ktoplot_1[i]) / Pnl_int(ktoplot_1[i]) / Pnl_int(squeezeness_1*ktoplot_1[i]))
    total_curve[i]        = 1.0e3*(       cov_total_mono_int(ktoplot_1[i],squeezeness_1*ktoplot_1[i]) / Pnl_int(ktoplot_1[i]) / Pnl_int(squeezeness_1*ktoplot_1[i]))
plot(ktoplot_1, stitchedtree_curve, linestyle = 'solid' , linewidth = 2.0, c = 'b', label = r"$Stitched$")
plot(ktoplot_1, oloop_curve       , linestyle = 'solid' , linewidth = 2.0, c = 'g', label = r"$1-loop$")
plot(ktoplot_1, total_curve       , linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"$Total$")

axvline(0.075/squeezeness_1, linestyle = 'dashed', c = 'k')
annotate(r"$II$"  ,xy=(0.55,0.75),xycoords='axes fraction',rotation=0,xytext=(0.55,0.75),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)
annotate(r"$III$" ,xy=(0.85,0.75),xycoords='axes fraction',rotation=0,xytext=(0.85,0.75),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)

xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ {\rm Cov}^{\rm NG, \ell=0}/(P_{m,1} P_{m,2})\ \times\ 10^{3}$'        , fontsize = labelsize)
#xlim(min((k_common)), max((k_common)))
xlim(0.04, 0.7)
ylim(0.0, 0.75)
xscale('log')
xticks(size = ticksize)
yticks(size = ticksize)
#params = {'legend.fontsize': legendsize}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 1)

fig0.add_subplot(2,2,2)

stitchedtree_curve = zeros(len(ktoplot_2))
oloop_curve        = zeros(len(ktoplot_2))
total_curve        = zeros(len(ktoplot_2))
for i in range(len(ktoplot_2)):
    stitchedtree_curve[i] = 1.0e3*(cov_stitchedtree_mono_int(ktoplot_2[i],squeezeness_2*ktoplot_2[i]) / Pnl_int(ktoplot_2[i]) / Pnl_int(squeezeness_2*ktoplot_2[i]))
    oloop_curve[i]        = 1.0e3*(       cov_1loop_mono_int(ktoplot_2[i],squeezeness_2*ktoplot_2[i]) / Pnl_int(ktoplot_2[i]) / Pnl_int(squeezeness_2*ktoplot_2[i]))
    total_curve[i]        = 1.0e3*(       cov_total_mono_int(ktoplot_2[i],squeezeness_2*ktoplot_2[i]) / Pnl_int(ktoplot_2[i]) / Pnl_int(squeezeness_2*ktoplot_2[i]))
plot(ktoplot_2, stitchedtree_curve, linestyle = 'solid' , linewidth = 2.0, c = 'b', label = r"$Stitched\ tree$")
plot(ktoplot_2, oloop_curve       , linestyle = 'solid' , linewidth = 2.0, c = 'g', label = r"$Response\ 1-loop$")
plot(ktoplot_2, total_curve       , linestyle = 'solid' , linewidth = 2.0, c = 'r', label = r"$Total$")

axvline(0.075/squeezeness_2, linestyle = 'dashed', c = 'k')
axvline(0.300/squeezeness_2, linestyle = 'dashed', c = 'k')
annotate(r"$I$"  ,xy=(0.20,0.45),xycoords='axes fraction',rotation=0,xytext=(0.20,0.45),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)
annotate(r"$IV$" ,xy=(0.55,0.45),xycoords='axes fraction',rotation=0,xytext=(0.55,0.45),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)
annotate(r"$V$"  ,xy=(0.90,0.45),xycoords='axes fraction',rotation=0,xytext=(0.90,0.45),textcoords='axes fraction',bbox=dict(boxstyle="round4",fc='w',alpha=0.75),color='k',fontsize = textsize)

xlabel(r'$k_1\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
#ylabel(r'${\rm Cov}_{12}/P_{m,1}/P_{m,2}\ \times\ 10^{3}$'        , fontsize = labelsize)
#xlim(min((k_common)), max((k_common)))
xlim(0.04, 0.7)
ylim(0.0, 0.75)
xscale('log')
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legendsize}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 2)

fig0.savefig('fig_regimes_v2.png')

show()
