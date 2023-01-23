#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 09:52:55 2022

@author: cmosbeux

Description: This script allows to make a visual comparison between different versions of BedMachine. 
The script could be adapted for a different comparison

The Bedmachine_version_difference.py uses classic python libraries/modules:
    - numpy
    - matplotlib
    - os
    
As well as other personal libraries available in the "Modules" subdirectory.
    
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


from Modules.format_reading import netcdf
from Modules.Antarctica_Background import Plot_Antarctica, basin

#%% definition of some functions

def create_directory(name="./Figures"):# checking if the directory demo_folder 
# exist or not.
    if not os.path.exists(name):
          
        # if the demo_folder directory is not present 
        # then create it.
        os.makedirs(name)


def symbolic_link(src):
    try:
        dst = './Data_SymLinks'
        create_directory(dst)
        dst_file = src.split('/')[-1]
        if not os.path.exists(src):
            print('Source file does not exist')
        os.symlink(src, dst+'/'+dst_file) 
        print("symbolic link to %s created in %s" % (dst_file, dst))
    except FileExistsError:
        message1 = "The symbolic link to %s already exists." % (dst_file)
        message2 = "Please, delete the link if you want to make a new one."
        print(message1, '\n', message2)
    
    
def slope(bed):
    slopes = np.gradient(bed)
    slope = (slopes[0] ** 2 + slopes[1] ** 2) ** 0.5
    return slope


def read_bedmachine(version="v1", resolution=1):
    print("Reading bedmachine...")
    print("\t - using a %0.1f km resolution" % (resolution / 2.0))

    process = True
    if version == "v1":
        bedmachine_filename = "./Data_SymLinks/BedMachineAntarctica_2019-11-05_v01.nc"
    elif version == "v2":
        bedmachine_filename = "./Data_SymLinks/BedMachineAntarctica_2020-07-15_v02.nc"
    elif version == "v3":
        bedmachine_filename = "./Data_SymLinks/BedMachineAntarctica-v03.nc"
    else:
        print("%s version is not available. Check that it does actually exist!" % version)
        process = False

    if process:
        x, y, bed = netcdf.readgrid(bedmachine_filename, "bed")
        x, y, surf = netcdf.readgrid(bedmachine_filename, "surface")
        x, y, thickness = netcdf.readgrid(bedmachine_filename, "thickness")
        x, y, mask = netcdf.readgrid(bedmachine_filename, "mask")

        x_inf, x_sup = +0, -1
        y_inf, y_sup = +0, -1
        x = x[x_inf:x_sup:resolution]
        y = y[y_inf:y_sup:resolution]
        bed = bed[x_inf:x_sup:resolution, y_inf:y_sup:resolution]
        surf = surf[x_inf:x_sup:resolution, y_inf:y_sup:resolution]
        thickness = thickness[x_inf:x_sup:resolution, y_inf:y_sup:resolution]
        mask = mask[x_inf:x_sup:resolution, y_inf:y_sup:resolution]

        mask2 = np.zeros_like(mask)
        mask2[mask == 0] = True
        mask2[mask > 0.5] = False

        bed = np.ma.masked_array(bed, mask2)
        surf = np.ma.masked_array(surf, mask2)
        thickness = np.ma.masked_array(thickness, mask2)

    return x, y, bed, surf, thickness


#%%
''' This section needs to be filled for the script to find the BedMachine Data'''

create_directory() #creation of a directory for the output figures
symbolic_link("/Users/cmosbeux/Documents/Data_Antarctica.CLEAN/DEMs&dHdT/BedMachineAntarctica_2019-11-05_v01.nc")
symbolic_link("/Users/cmosbeux/Documents/Data_Antarctica.CLEAN/DEMs&dHdT/BedMachineAntarctica_2020-07-15_v02.nc")
symbolic_link("/Users/cmosbeux/Documents/Data_Antarctica.CLEAN/DEMs&dHdT/BedMachineAntarctica-v03.nc")


'''The resolution of the grid can be changed to speed up the process. WARNING: resolution=1 (i.e. 0.5 km) is for native resolution.
A higher number will deteriorate the evalution and should only be used for quick tests'''

resolution=5

#%% Start grid comparisons and grid calculations 
x, y, bed_v1, surf, thickness_v1 = read_bedmachine(version="v1", resolution=resolution)
x, y, bed_v2, surf, thickness_v2 = read_bedmachine(version="v2", resolution=resolution)
x, y, bed_v3, surf, thickness_v3 = read_bedmachine(version="v3", resolution=resolution)

#get the gradients
grad_v1 = slope(bed_v1)
grad_v2 = slope(bed_v2)
grad_v3 = slope(bed_v3)


#get the map extent
extent_IS = [np.min(x), np.max(x), np.min(y), np.max(y)]
    

#%%

params = {
    "ytick.color": "k",
    "xtick.color": "k",
    "axes.labelcolor": "k",
    "axes.edgecolor": "k",
}
plt.rcParams.update(params)

plot = True
plot_panantarctic = True


if plot:
    extent = np.asarray(basin.Amundsen())
    fig, ax = Plot_Antarctica(
        nrows=1,
        ncols=1,
        basemap="light",
        continental_shelf=0.0,
        extent=extent,
        figsize=(30, 25),
    )

    cb = ax[0].imshow(
        thickness_v3 - thickness_v2,
        cmap="seismic",
        extent=extent_IS,
        vmin=-200,
        vmax=200,
        alpha=1,
        zorder=1e1,
    )

    cbar = ax[0].cax.colorbar(cb)
    cbar = ax.cbar_axes[0].colorbar(cb)
    cbaxes = fig.add_axes([0.33, 0.175, 0.12, 0.015])
    cbar = fig.colorbar(cb, cax=cbaxes, ticks=[-200, 0, 200], orientation="horizontal")
    cbar.set_label("ice thickness difference", fontsize=16)
    cbar.ax.set_xticklabels(["-200", "0", r"   200 m"], fontsize=13)
    plt.savefig(
        "Bedmachine_Thickness_V3-V2.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.1,
        transparent=False,
    )

if plot:
    extent = np.asarray(basin.Amundsen())
    fig, ax = Plot_Antarctica(
        nrows=1,
        ncols=1,
        basemap="light",
        continental_shelf=0.0,
        extent=extent,
        figsize=(30, 25),
    )

    cb = ax[0].imshow(
        bed_v3 - bed_v2,
        cmap="seismic",
        extent=extent_IS,
        vmin=-200,
        vmax=200,
        alpha=1,
        zorder=1e1,
    )

    cbar = ax[0].cax.colorbar(cb)
    cbar = ax.cbar_axes[0].colorbar(cb)
    cbaxes = fig.add_axes([0.33, 0.175, 0.12, 0.015])
    cbar = fig.colorbar(cb, cax=cbaxes, ticks=[-200, 0, 200], orientation="horizontal")
    cbar.set_label("Bed difference", fontsize=16)
    cbar.ax.set_xticklabels(["-200", "0", r"   200 m"], fontsize=13)
    plt.savefig(
        "Figures/Bedmachine_Bed_V3-V2.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.1,
        transparent=False,
    )

plot_grad = True
if plot_grad:
    extent = np.asarray(basin.Amundsen())
    fig, ax = Plot_Antarctica(
        nrows=1,
        ncols=1,
        basemap="light",
        continental_shelf=0.0,
        extent=extent,
        figsize=(30, 25),
    )

    cb = ax[0].imshow(
        grad_v2 * 1e-3,
        cmap="hot",
        extent=extent_IS,
        norm=mpl.colors.LogNorm(vmin=0.01, vmax=1),
        alpha=1,
        zorder=1e1,
    )

    cbar = ax[0].cax.colorbar(cb)
    cbar = ax.cbar_axes[0].colorbar(cb)
    cbaxes = fig.add_axes([0.33, 0.175, 0.12, 0.015])
    cbar = fig.colorbar(cb, cax=cbaxes, ticks=[1e-2, 1e-1, 1], orientation="horizontal")
    cbar.set_label("Bed Slope", fontsize=16)
    cbar.ax.set_xticklabels([r"10$^{-2}$", r"10$^{-1}$", r"1"], fontsize=13)
    plt.savefig(
        "Figures/Bedmachine_BedSlope_v2.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.1,
        transparent=False,
    )

    fig, ax = Plot_Antarctica(
        nrows=1,
        ncols=1,
        basemap="light",
        continental_shelf=0.0,
        extent=extent,
        figsize=(30, 25),
    )

    cb = ax[0].imshow(
        grad_v3 * 1e-3,
        cmap="hot",
        extent=extent_IS,
        norm=mpl.colors.LogNorm(vmin=0.01, vmax=1),
        alpha=1,
        zorder=1e1,
    )

    cbar = ax[0].cax.colorbar(cb)
    cbar = ax.cbar_axes[0].colorbar(cb)
    cbaxes = fig.add_axes([0.33, 0.175, 0.12, 0.015])
    cbar = fig.colorbar(cb, cax=cbaxes, ticks=[1e-2, 1e-1, 1], orientation="horizontal")
    cbar.set_label("Bed Slope", fontsize=16)
    cbar.ax.set_xticklabels([r"10$^{-2}$", r"10$^{-1}$", r"  1"], fontsize=13)
    plt.savefig(
        "Figures/Bedmachine_BedSlope_v3.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.1,
        transparent=False,
    )

    fig, ax = Plot_Antarctica(
        nrows=1,
        ncols=1,
        basemap="light",
        continental_shelf=0.0,
        extent=extent,
        figsize=(30, 25),
    )

    dgrad = (grad_v3 - grad_v2) * 1e-3
    cb = ax[0].imshow(
        dgrad,
        cmap="seismic",
        extent=extent_IS,
        norm=mpl.colors.Normalize(vmin=-1e-1, vmax=1e-1),
        alpha=1,
        zorder=1e1,
    )

    cbar = ax[0].cax.colorbar(cb)
    cbar = ax.cbar_axes[0].colorbar(cb)
    cbaxes = fig.add_axes([0.33, 0.175, 0.12, 0.015])
    cbar = fig.colorbar(
        cb, cax=cbaxes, ticks=[-1e-1, 0, 1e-1], orientation="horizontal"
    )
    cbar.set_label("Bed Slope difference_V3-V2", fontsize=16)
    cbar.ax.set_xticklabels([r"10$^{-1}$", r"0", r"10$^{-1}$"], fontsize=13)
    plt.savefig(
        "Figures/Bedmachine_BedSlope_v3-v2.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.1,
        transparent=False,
    )


if plot_panantarctic:
    extent = basin.PanAntarctic()
    fig, ax = Plot_Antarctica(
        nrows=1,
        ncols=1,
        basemap="light",
        continental_shelf=0.0,
        extent=extent,
        figsize=(30, 25),
    )

    cb = ax[0].imshow(
        thickness_v3 - thickness_v2,
        cmap="seismic",
        extent=extent_IS,
        vmin=-200,
        vmax=200,
        alpha=1,
        zorder=1e1,
    )
    cbar = ax[0].cax.colorbar(cb)
    cbar = ax.cbar_axes[0].colorbar(cb)
    cbaxes = fig.add_axes([0.35, 0.175, 0.10, 0.015])
    cbar = fig.colorbar(cb, cax=cbaxes, ticks=[-200, 0, 200], orientation="horizontal")
    cbar.set_label("ice thickness difference", fontsize=14)
    cbar.ax.set_xticklabels(["-200", "0", r"   200 m"], fontsize=13)
    plt.savefig(
        "Figures/Bedmachine_Thickness_PA_V3-V2.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.1,
        transparent=False,
    )
