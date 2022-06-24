#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 11:15:01 2021

@author: ltomkins
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from utils import gen_fun
import cmocean
import numpy as np
import pyart

def create_plotDict(var):
    """
    Creates a dictionary used in the plotting functions

    Parameters
    ----------
    var : str - Name of variable to return plotting information for

    Returns
    -------
    plot_dict : dict - dictionary of plotting information
    """
    magma_cmap = plt.get_cmap('magma_r')
    grays_cmap = plt.get_cmap('gray_r')
    rdbu_cmap = plt.get_cmap('RdBu_r')

    if var in ('ref_max', 'reflectivity', 'nonmuted_ref', 'nonmuted_reflectivity'):
        plot_dict = {}
        plot_dict['colormap'] = mcolors.LinearSegmentedColormap.from_list('nonmuted_cmap', magma_cmap(np.linspace(0, 0.9, magma_cmap.N)))
        plot_dict['vmin'] = 5; plot_dict['vmax'] = 45;
        plot_dict['title'] = 'Reflectivity [dBZ]'
        plot_dict['tick_flag'] = False
        plot_dict['tick_label_flag'] = False

    elif var in ('muted_ref', 'muted_ref_max', 'muted_reflectivity'):
        plot_dict = {}
        plot_dict['colormap'] = mcolors.LinearSegmentedColormap.from_list('muted_cmap', grays_cmap(np.linspace(0, 0.7, grays_cmap.N)))
        plot_dict['vmin'] = 5;  plot_dict['vmax'] = 45;
        plot_dict['title'] = 'Muted Reflectivity [dBZ]'
        plot_dict['tick_flag'] = False
        plot_dict['tick_label_flag'] = False

    elif var in ('waves_max'):
        plot_dict = {}
        plot_dict['colormap'] = mcolors.LinearSegmentedColormap.from_list('waves', [(1, 1, 1), (0, 0, 0)], N=2)
        plot_dict['vmin'] = 0; plot_dict['vmax'] = 1;
        plot_dict['title'] = 'Waves'
        plot_dict['tick_flag'] = True; plot_dict['ticks'] = [0.25, 0.75]
        plot_dict['tick_label_flag'] = True; plot_dict['tick_labels'] = ['', 'Wave']

    elif var in ('rho_max', 'cross_correlation_ratio'):
        plot_dict = {}
        plot_dict['colormap'] = cmocean.tools.crop_by_percent(cmocean.cm.deep, 10, which='max',N=None)
        plot_dict['vmin'] = 0.8; plot_dict['vmax'] = 1
        plot_dict['title'] = 'Correlation Coefficient'
        plot_dict['tick_flag'] = True; plot_dict['ticks'] = [0.8, 0.85, 0.9, 0.95, 0.97, 1.0]
        plot_dict['tick_label_flag'] = False

    elif var in ('velocity'):
        plot_dict = {}
        plot_dict['colormap'] = rdbu_cmap
        plot_dict['vmin'] = -40; plot_dict['vmax'] = 40
        plot_dict['title'] = 'Velocity [ms$^{-1}$]'
        plot_dict['tick_flag'] = False
        plot_dict['tick_label_flag'] = False

    else:
        print(var + ' Key not recognized')

    return plot_dict

def xytoR(x, y, azi_mean):
    """
    Calculates R used to plot RHIs

    Parameters
    ----------
    x : array - Array of x values
    y : array - Array of y values
    azi_mean : float - value of average azimuth

    Returns
    -------
    R : array - array of R values
    """

    if 89.5 <= azi_mean <= 90.0:
        R = np.sqrt(x ** 2 + y ** 2) * np.sign(x)
    else:
        R = np.sqrt(x ** 2 + y ** 2) * np.sign(y)
    return R

def plotRHI(ax, azimuth, display_rhi, plot_var, muted=False, muted_var = None,
            aspect_ratio=5, xlims=(0,200), ylims=(0,6),title_meta_flag=True,
            title_var_flag=True, colorbar_flag=True, xax_label=True, yax_label=True):
    """
    This function plots RHIs
    :param ax: axes object to plot on
    :param azimuth: value of RHI azimuth
    :param display_rhi: pyart display object with RHI information
    :param plot_var: str - variable to plot
    :param muted: bool - to mute (True) variable or not (False)
    :param muted_var: str - muted variable
    :param aspect_ratio: int - aspect ratio
    :param xlims: tuple of ints - x axis limits
    :param ylims: tuples of ints - y axis limits
    :param title_meta_flag: bool - to show metadata tile
    :param title_var_flag: bool - to show variable title
    :param colorbar_flag: bool - to plot colorbar
    :param xax_label: bool - to plot xaxis labels
    :param yax_label: bool - to plot yaxis labels
    :return: ax with plot
    """

    plot_dict = create_plotDict(plot_var)

    data = display_rhi._get_data(plot_var, 0, None, False, None)
    x, y, z = display_rhi._get_x_y_z(0, False, False)
    azi_mean = np.abs(np.mean(display_rhi._radar.azimuth['data'][display_rhi._radar.get_slice(0)]))
    R = xytoR(x, y, azi_mean)

    # plot main variable
    pm = ax.pcolormesh(R, z, data, cmap=plot_dict['colormap'], vmin=plot_dict['vmin'],
                       vmax=plot_dict['vmax'], shading='nearest', rasterized=True)

    # set colorbar
    if colorbar_flag:
        cbar_pad_val = 0.01
        if muted: cbar_pad_val = -0.1
        if not plot_dict['tick_flag']: cbar = plt.colorbar(pm, ax=ax, fraction=0.040, pad=cbar_pad_val)
        if plot_dict['tick_flag']: cbar = plt.colorbar(pm, ax=ax, fraction=0.040, pad=cbar_pad_val, ticks=plot_dict['ticks'])
        if plot_dict['tick_label_flag']: cbar.ax.set_yticklabels(plot_dict['tick_labels'])

    # plot muted variable
    if muted:
        muted_dict = create_plotDict(muted_var)
        muted_data = display_rhi._get_data(muted_var, 0, None, False, None)
        pm_mask = ax.pcolormesh(R, z, muted_data, cmap=muted_dict['colormap'], vmin=muted_dict['vmin'],
                            vmax=muted_dict['vmax'], shading='nearest', rasterized=True)
        if colorbar_flag:
            cbar1 = plt.colorbar(pm_mask, ax=[ax], pad=0.01);
            cbar1.ax.tick_params(labelsize=0)

    # set aspect ratio
    ax.set_aspect(aspect_ratio)

    # set plot limits
    if np.mean(R) <= 0: xlims = (xlims[0]*-1, xlims[1]*-1)
    ax.set(xlim=xlims, ylim=ylims)
    xvals = np.arange(np.floor(xlims[0] / 25), np.floor(xlims[1] / 25)*25 + 25, 25)
    yvals = np.arange(np.floor(ylims[0]), np.floor(ylims[1]/2)*2+2, 2)
    ax.set_xticks(xvals)
    ax.set_yticks(yvals)

    # set axis labels
    if xax_label: ax.set_ylabel('Distance from radar [km]')
    if yax_label: ax.set_ylabel('Height [km]')

    # set titles
    if title_meta_flag:
        ax.set_title(u'{0} {1}\N{DEGREE SIGN} azimuth cross-section \n{2}'.
                     format(pyart.graph.common.generate_radar_name(display_rhi._radar), azimuth,
                            pyart.util.datetime_from_radar(display_rhi._radar).strftime('%d %b %Y %H:%M:%S UTC')), loc='left')

    if title_var_flag: ax.set_title(plot_dict['title'], loc='right')

    return ax


def plotPPI(ax, display_ppi, plot_var, elev_angle=0.5, muted=False, muted_var = None,
            xlims=(-200,200), ylims=(-200, 200),title_flag=False,
            colorbar_flag=True, xax_label=True, yax_label=True,
            plot_RHI_loc=False, azimuth=None, RHI_loc_color=None):
    """
    This function plots PPIs
    :param ax: axes object to plot on
    :param display_ppi: pyart display object with PPI information
    :param plot_var: str - variable to plot
    :param elev_angle: float - elevation angle to plot
    :param muted: bool - to mute (True) variable or not (False)
    :param muted_var: str - muted variable
    :param xlims: tuple of ints - x axis limits
    :param ylims: tuples of ints - y axis limits
    :param title_flag: bool - to show metadata tile
    :param colorbar_flag: bool - to plot colorbar
    :param xax_label: bool - to plot xaxis labels
    :param yax_label: bool - to plot yaxis labels
    :param plot_RHI_loc: bool - to plot location of an RHI
    :param azimuth: value of RHI azimuth
    :param RHI_loc_color: str - color of azimuth
    :return: ax with plot
    """

    plot_dict = create_plotDict(plot_var)

    data = display_ppi._get_data(plot_var, 0, None, False, None)
    x, y, _ = display_ppi._get_x_y_z(0, False, False)

    # plot main variable
    pm = ax.pcolormesh(x, y, data, cmap=plot_dict['colormap'], vmin=plot_dict['vmin'],
                       vmax=plot_dict['vmax'], shading='nearest', rasterized=True)

    # set colorbar
    if colorbar_flag:
        cbar_pad_val = 0.01
        if muted: cbar_pad_val = -0.1
        if not plot_dict['tick_flag']: cbar = plt.colorbar(pm, ax=ax, fraction=0.040, pad=cbar_pad_val)
        if plot_dict['tick_flag']: cbar = plt.colorbar(pm, ax=ax, fraction=0.040, pad=cbar_pad_val, ticks=plot_dict['ticks'])
        if plot_dict['tick_label_flag']: cbar.ax.set_yticklabels(plot_dict['tick_labels'])

    # plot muted variable
    if muted:
        muted_dict = create_plotDict(muted_var)
        muted_data = display_ppi._get_data(muted_var, 0, None, False, None)
        pm_mask = ax.pcolormesh(x, y, muted_data, cmap=muted_dict['colormap'], vmin=muted_dict['vmin'],
                            vmax=muted_dict['vmax'], shading='nearest', rasterized=True)
        if colorbar_flag:
            cbar1 = plt.colorbar(pm_mask, ax=[ax], pad=0.01);
            cbar1.ax.tick_params(labelsize=0)

    # plot RHI location
    if plot_RHI_loc:
        x_c, y_c = gen_fun.azi_calculator([azimuth], np.hypot(xlims[1], ylims[1]))
        ax.plot([0, x_c[0]], [0, y_c[0]], color=RHI_loc_color)

    # set aspect ratio
    ax.set_aspect('equal')

    # set plot limits
    display_ppi.set_limits(xlim=xlims, ylim=ylims)
    xvals = np.arange(np.floor(xlims[0] / 100)*100, np.floor(xlims[1] / 100)*100 + 100, 100)
    yvals = np.arange(np.floor(ylims[0] / 100)*100, np.floor(ylims[1] / 100)*100 + 100, 100)
    ax.set_xticks(xvals)
    ax.set_yticks(yvals)

    # set axis labels
    if xax_label: ax.set_ylabel('Distance from radar [km]')
    if yax_label: ax.set_ylabel('Distance from radar [km]')

    # set titles
    if title_flag:
        ax.set_title(u'{0} {1}\N{DEGREE SIGN} tilt {2} \n{3}'.
                     format(pyart.graph.common.generate_radar_name(display_ppi._radar), elev_angle, plot_dict['title'],
                            pyart.util.datetime_from_radar(display_ppi._radar).strftime('%d %b %Y %H:%M:%S UTC')), loc='left')


    return ax


# scratch
#
#     elif field == 'snow_rate':
#         sn_dict = {}
#         sn_dict['data'] = grid.fields[field]['data'][0]
#         sn_dict['colormap'] = 'viridis'
#         sn_dict['vmin'] = 0;
#         sn_dict['vmax'] = 5;
#         sn_dict['title'] = 'Snow Rate [mm hr$^{-1}$]'
#         sn_dict['tick_flag'] = False
#         sn_dict['tick_label_flag'] = False
#
#         plot_dict[field] = sn_dict
#
#     elif field in ('convsf', 'convsf2'):
#         convsf_dict = {}
#         convsf_dict['data'] = grid.fields[field]['data'][0]
#         convsf_dict['colormap'] = LinearSegmentedColormap.from_list('convsf', [(1, 1, 1), (0.1287, 0.5633, 0.5512),
#                                                                                (0.9932, 0.9062, 0.1439)], N=3)
#         convsf_dict['vmin'] = 0;
#         convsf_dict['vmax'] = 2;
#         convsf_dict['title'] = ''
#         convsf_dict['tick_flag'] = True;
#         convsf_dict['ticks'] = [0.25, 1, 1.75]
#         convsf_dict['tick_label_flag'] = True;
#         convsf_dict['tick_labels'] = ['', 'Stratiform', 'Convective']
#
#         plot_dict[field] = convsf_dict
#
#     elif field == 'alt_max':
#         alt_dict = {}
#         alt_dict['data'] = grid.fields[field]['data'][0]
#         alt_dict['colormap'] = 'cmo.haline'
#         alt_dict['vmin'] = 0;
#         alt_dict['vmax'] = 5000;
#         alt_dict['title'] = 'Altitude [m]'
#         alt_dict['tick_flag'] = False
#         alt_dict['tick_label_flag'] = False
#
#         plot_dict[field] = alt_dict
