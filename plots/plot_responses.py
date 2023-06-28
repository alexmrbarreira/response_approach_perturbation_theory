import sys; sys.path.append('../'); 
from parameters import *
from functions import *

# ==================================================================== #
# Plot the response functions 
# ==================================================================== #

font = 18
fontlabel = 22
legend_font = 18
text_font = 20
marker_size = 40
marker_edge_width = 1.
line_width = 1.5
x_text = 2.0e14
y_text = 10.3
x_min = 5.0e-3
#x_max = 6.
x_max = 30
tick_length = 4.0
tick_width  = 2.0
binsy = 6
c_1  = 'b'
c_2  = 'g'
c_3  = 'r'
c_4  = 'k'
c_5  = 'b'
c_6  = 'g'
c_7  = 'r'
c_8  = 'k'
ls_1 = 'solid'
ls_2 = 'solid'
ls_3 = 'solid'
ls_4 = 'solid'
ls_5 = 'dashed'
ls_6 = 'dashed'
ls_7 = 'dashed'
ls_8 = 'dashed'
lab_1  = r"$R_1(k)$"
lab_2  = r"$R_K(k)$"
lab_3  = r"$R_2(k)/2$"
lab_4  = r"$R_{K\delta}(k)$"
lab_5  = r"$R_{K^2}(k)$"
lab_6  = r"$R_{K.K}(k)$"
lab_7  = r"$R_{KK}(k)$"
lab_8  = r"$R_{\hat{\Pi}}(k)$"
lab_1t = r"$R_1^{L,\rm tree}(k)$"
lab_2t = r"$R_K^{L,\rm tree}(k)$"
lab_3t = r"$R_2^{L,\rm tree}(k)$"
lab_4t = r"$R_{K\delta}^{L,\rm tree}(k)$"
lab_5t = r"$R_{K^2}^{L,\rm tree}(k)$"
lab_6t = r"$R_{K.K}^{L,\rm tree}(k)$"
lab_7t = r"$R_{KK}^{L,\rm tree}(k)$"
lab_8t = r"$R_{\hat{\Pi}}^{L,\rm tree}(k)$"
lab_1n = r"$R_1^{L,\rm nl}(k)$"
lab_2n = r"$R_K^{L,\rm nl}(k)$"
lab_3n = r"$R_2^{L,\rm nl}(k)$"
lab_4n = r"$R_{K\delta}^{L,\rm nl}(k)$"
lab_5n = r"$R_{K^2}^{L,\rm nl}(k)$"
lab_6n = r"$R_{K.K}^{L,\rm nl}(k)$"
lab_7n = r"$R_{KK}^{L,\rm nl}(k)$"
lab_8n = r"$R_{\hat{\Pi}}^{L,\rm nl}(k)$"

# ==================================================================== #
# Eulerian figure
# ==================================================================== #
fig0 = plt.figure(0, figsize=(17., 6.))
fig0.subplots_adjust(left=0.06, right=0.98, top=0.96, bottom= 0.14)#, hspace = 0.15)#, space = 0.2)
fig0.add_subplot(1,2,1)
plot(k_common, R1e_tree(k_common)   , c = c_1, linestyle = ls_1, linewidth = line_width, label = lab_1)
plot(k_common, R2e_tree(k_common)   , c = c_2, linestyle = ls_2, linewidth = line_width, label = lab_2)
plot(k_common, R3e_tree(k_common)/2., c = c_3, linestyle = ls_3, linewidth = line_width, label = lab_3)
plot(k_common, R4e_tree(k_common)   , c = c_4, linestyle = ls_4, linewidth = line_width, label = lab_4)
plot(k_common, R5e_tree(k_common)   , c = c_5, linestyle = ls_5, linewidth = line_width, label = lab_5)
plot(k_common, R6e_tree(k_common)   , c = c_6, linestyle = ls_6, linewidth = line_width, label = lab_6)
plot(k_common, R7e_tree(k_common)   , c = c_7, linestyle = ls_7, linewidth = line_width, label = lab_7)
plot(k_common, R8e_tree(k_common)   , c = c_8, linestyle = ls_8, linewidth = line_width, label = lab_8)
annotate(r"$Tree\ level$", xy = (0.04,0.57), xycoords='axes fraction', color = 'k', fontsize = text_font+1)
xscale('log')
xlabel(r"$k\ \left[h/{\rm Mpc}\right]$", fontsize = fontlabel)
ylabel(r"$Response\ coefficients,\ R_\mathcal{O}(k)$", fontsize = fontlabel)
xlim(min(k_common), max(k_common))
ylim(ymin = -7, ymax = 25)
xticks(size = font)
yticks(size = font)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 3)
fig0.add_subplot(1,2,2)
plot(k_common, R1e_nonl(k_common)   , c = c_1, linestyle = ls_1, linewidth = line_width, label = lab_1)
plot(k_common, R2e_nonl(k_common)   , c = c_2, linestyle = ls_2, linewidth = line_width, label = lab_2)
plot(k_common, R3e_nonl(k_common)/2., c = c_3, linestyle = ls_3, linewidth = line_width, label = lab_3)
plot(k_common, R4e_nonl(k_common)   , c = c_4, linestyle = ls_4, linewidth = line_width, label = lab_4)
plot(k_common, R5e_nonl(k_common)   , c = c_5, linestyle = ls_5, linewidth = line_width, label = lab_5)
plot(k_common, R6e_nonl(k_common)   , c = c_6, linestyle = ls_6, linewidth = line_width, label = lab_6)
plot(k_common, R7e_nonl(k_common)   , c = c_7, linestyle = ls_7, linewidth = line_width, label = lab_7)
plot(k_common, R8e_nonl(k_common)   , c = c_8, linestyle = ls_8, linewidth = line_width, label = lab_8)
annotate(r"$Nonlinear$", xy = (0.04,0.57), xycoords='axes fraction', color = 'k', fontsize = text_font+1)
xscale('log')
xlabel(r"$k\ \left[h/{\rm Mpc}\right]$", fontsize = fontlabel)
ylabel(r"$Response\ coefficients,\ R_\mathcal{O}(k)$", fontsize = fontlabel)
xlim(min(k_common), max(k_common))
ylim(ymin = -7, ymax = 25)
xticks(size = font)
yticks(size = font)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 3)

# ==================================================================== #
# Eulerian vs Lagrangian figure
# ==================================================================== #
fig2 = plt.figure(2, figsize=(17., 6.))
fig2.subplots_adjust(left=0.06, right=0.98, top=0.96, bottom= 0.14)#, hspace = 0.15)#, space = 0.2)
fig2.add_subplot(1,2,1)
plot(k_common, R3l_tree(k_common), c = 'b', linestyle = 'solid', linewidth = line_width, label = lab_3)
plot(k_common, R4l_tree(k_common), c = 'g', linestyle = 'solid', linewidth = line_width, label = lab_4)
plot(k_common, R5l_tree(k_common), c = 'r', linestyle = 'solid', linewidth = line_width, label = lab_5)
plot(k_common, R6l_tree(k_common), c = 'c', linestyle = 'solid', linewidth = line_width, label = lab_6)
plot(k_common, R8l_tree(k_common), c = 'k', linestyle = 'solid', linewidth = line_width, label = lab_8)
plot(k_common, R3e_tree(k_common), c = 'b', linestyle = 'dashed', linewidth = line_width)
plot(k_common, R4e_tree(k_common), c = 'g', linestyle = 'dashed', linewidth = line_width)
plot(k_common, R5e_tree(k_common), c = 'r', linestyle = 'dashed', linewidth = line_width)
plot(k_common, R6e_tree(k_common), c = 'c', linestyle = 'dashed', linewidth = line_width)
plot(k_common, R8e_tree(k_common), c = 'k', linestyle = 'dashed', linewidth = line_width)
annotate(r"$Tree\ level$"         , xy = (0.02,0.68)         , xycoords='axes fraction', color = 'k', fontsize = text_font+1)
annotate(r"$Solid\ \ \ \ \ :\ Lagrangian$", xy = (0.02,0.59), xycoords='axes fraction', color = 'k', fontsize = text_font-1)
annotate(r"$Dashed\ :\ Eulerian$" , xy = (0.02,0.52) , xycoords='axes fraction', color = 'k', fontsize = text_font-1)
xscale('log')
xlabel(r"$k\ \left[h/{\rm Mpc}\right]$", fontsize = fontlabel)
ylabel(r"$Response\ functions$", fontsize = fontlabel)
xlim(min(k_common), max(k_common))
ylim(ymin = -7, ymax = 25)
xticks(size = font)
yticks(size = font)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 3)
fig2.add_subplot(1,2,2)
plot(k_common, R3l_nonl(k_common), c = 'b', linestyle = 'solid', linewidth = line_width, label = lab_3)
plot(k_common, R4l_nonl(k_common), c = 'g', linestyle = 'solid', linewidth = line_width, label = lab_4)
plot(k_common, R5l_nonl(k_common), c = 'r', linestyle = 'solid', linewidth = line_width, label = lab_5)
plot(k_common, R6l_nonl(k_common), c = 'c', linestyle = 'solid', linewidth = line_width, label = lab_6)
plot(k_common, R8l_nonl(k_common), c = 'k', linestyle = 'solid', linewidth = line_width, label = lab_8)
plot(k_common, R3e_nonl(k_common), c = 'b', linestyle = 'dashed', linewidth = line_width)
plot(k_common, R4e_nonl(k_common), c = 'g', linestyle = 'dashed', linewidth = line_width)
plot(k_common, R5e_nonl(k_common), c = 'r', linestyle = 'dashed', linewidth = line_width)
plot(k_common, R6e_nonl(k_common), c = 'c', linestyle = 'dashed', linewidth = line_width)
plot(k_common, R8e_nonl(k_common), c = 'k', linestyle = 'dashed', linewidth = line_width)
annotate(r"$Nonlinear$"         , xy = (0.02,0.68)         , xycoords='axes fraction', color = 'k', fontsize = text_font+1)
annotate(r"$Solid\ \ \ \ \ :\ Lagrangian$", xy = (0.02,0.59), xycoords='axes fraction', color = 'k', fontsize = text_font-1)
annotate(r"$Dashed\ :\ Eulerian$" , xy = (0.02,0.52) , xycoords='axes fraction', color = 'k', fontsize = text_font-1)
xscale('log')
xlabel(r"$k\ \left[h/{\rm Mpc}\right]$", fontsize = fontlabel)
ylabel(r"$Response\ functions$", fontsize = fontlabel)
xlim(min(k_common), max(k_common))
ylim(ymin = -7, ymax = 25)
xticks(size = font)
yticks(size = font)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 3)

fig0.savefig('fig_responses_eul.png')
fig2.savefig('fig_responses_eulvslag.png')

show()
