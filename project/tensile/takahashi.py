import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import os
import shutil
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from matplotlib.colors import Normalize
from matplotlib.colors import ListedColormap
from scipy.interpolate import griddata
from tqdm import tqdm

plt.rcParams["font.family"] ="DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "stix"

outputDirPath           = "/home/wtakahashi/plot/figure"
outputDirName           = "colloquium"
outputFileName          = "Kernel_form_1D_eta_1_02.png"
dir_path                = "/home/wtakahashi/result/result_boussinesq/"

filename_Wendland_W_1     = "boussinesq04_Kernel_check_params_form_Wendland_04.dat/boussinesq04_Wendland_1D_04_params_form_Wendland_04.dat_W.dat"
filename_Gaussian_W_1     = "boussinesq04_Kernel_check_params_form_Gaussian_04.dat/boussinesq04_Gaussian_1D_04_params_form_Gaussian_04.dat_W.dat"
filename_Spline_W_1       = "boussinesq04_Kernel_check_params_form_Spline_04.dat/boussinesq04_Spline_1D_04_params_form_Spline_04.dat_W.dat"

def plotData(dir_path_str, filename_str, x_label, y_label, data_label, data_color, data_linestyle, data_linewidth, ax):
    filename = os.path.join(dir_path_str, filename_str)
    with open(filename, "r") as file:
        data = np.loadtxt(file)
    x_arg = data[:, 0]
    y_arg = data[:, 1]
    ax.plot(x_arg, y_arg, label=data_label, color=data_color, linestyle=data_linestyle, linewidth=data_linewidth)
    ax.tick_params(labelsize=18)
    ax.tick_params(labelsize=18)
    ax.set_xlabel(x_label, fontsize=22)
    ax.set_ylabel(y_label, fontsize=22)

fig = plt.figure(figsize=(9, 6), tight_layout=True)
ax  = fig.add_subplot(111)

plotData(dir_path, filename_Wendland_W_1    , r"$x / h$", r"$W(x / h)$", r"$\mathrm{Wendland\ C4}$"  , "r", "-" , 2, ax)
plotData(dir_path, filename_Spline_W_1      , r"$x / h$", r"$W(x / h)$", r"$\mathrm{Cubic\ Spline}$" , "b", "-" , 2, ax)
plotData(dir_path, filename_Gaussian_W_1    , r"$x / h$", r"$W(x / h)$", r"$\mathrm{Gaussian}$"     , "g", "-" , 2, ax)

ax.set_xlim([-1.3,1.3])
ax.set_xlim([-3.0,3.0])
ax.set_ylim([0.0,1.55])
ax.legend(fontsize=18)

outputFilePath = os.path.join(outputDirPath, outputDirName, outputFileName)
fig.savefig(outputFilePath, transparent=True)
# plt.show()
