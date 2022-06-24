#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 11:11:55 2021

@author: ltomkins
"""

from utils import rhi_fun
from utils import gen_fun
from utils import plotRHI_fun
import pyart
import matplotlib.pyplot as plt
from datetime import datetime

#%% Variables to adjust
start_datetime = datetime.strptime('20191201_170000', '%Y%m%d_%H%M%S')
end_datetime   = datetime.strptime('20191201_200000', '%Y%m%d_%H%M%S')
radar          = 'KDIX'
azimuth        = 315

#%%
# file list
files = gen_fun.get_AWSlist(start_datetime, end_datetime, radar)

# read in radar object (PPI) - first file
radar = pyart.io.read_nexrad_archive(files[0])

# mute radar objects
radar = pyart.util.image_mute_radar(radar, 'reflectivity', 'cross_correlation_ratio', 0.97, field_threshold=20)

# RHI object
xsect = rhi_fun.constructRHI(radar, azimuth, 'reflectivity')

# Plot
display_ppi = pyart.graph.RadarDisplay(radar)
display_rhi = pyart.graph.RadarDisplay(xsect)

#%% using new plotting functions

# plot PPI
fig = plt.figure()
ax = plt.axes()
plotRHI_fun.plotPPI(ax, azimuth, display_ppi, 'cross_correlation_ratio', elev_angle=0.5, muted=False, muted_var = 'muted_reflectivity',
            xlims=(-200,200), ylims=(-200, 200),title_flag=True,
            colorbar_flag=True, xax_label=True, yax_label=True,
            plot_RHI_loc=True, RHI_loc_color='limegreen')
plt.show()

# plot RHI
fig = plt.figure(figsize=(10, 2))
ax = plt.axes()
plotRHI_fun.plotRHI(ax, azimuth, display_rhi, 'cross_correlation_ratio', muted=False, muted_var = 'muted_reflectivity',
            aspect_ratio=5, xlims=(0,200), ylims=(0,6),title_meta_flag=True,
            title_var_flag=True, colorbar_flag=True, xax_label=True, yax_label=True)
plt.show()


#%%
# adjust colormaps for visual separation
# this example uses perceptually uniform colormaps
# magma_cmap = plt.get_cmap('magma_r')
# grays_cmap = plt.get_cmap('gray_r')
#
# nonmuted_cmap = mcolors.LinearSegmentedColormap.from_list('nonmuted_cmap', magma_cmap(np.linspace(0, 0.9, magma_cmap.N)))
# muted_cmap = mcolors.LinearSegmentedColormap.from_list('muted_cmap', grays_cmap(np.linspace(0, 0.7, grays_cmap.N)))
#
# fig = plt.figure(figsize=(10, 3))
# ax = plt.axes()
#
# nm_ref = display_rhi.plot('nonmuted_reflectivity', vmin=5, vmax=45, colorbar_flag=False, title_flag=False,
#                  cmap=nonmuted_cmap)#, colorbar_orient='horizontal',colorbar_label='', ticks=[5,15,25,35,45])
# m_ref = display_rhi.plot('muted_reflectivity', vmin=5, vmax=45, colorbar_flag=False, title_flag=False,
#                  cmap=muted_cmap)#, colorbar_orient='horizontal',colorbar_label='', ticks=[5,15,25,35,45], ticklabs=None)
# display_rhi.set_aspect_ratio(5)
# display_rhi.set_limits(xlim=(0, 200), ylim=(0, 6))
# #cbar2 = plt.colorbar(nm_ref, ax=ax, pad=-0.08, ticks=[5,15,25,35,45])
# #cbar1 = plt.colorbar(m_ref, ax=ax, pad=0.02); cbar1.ax.tick_params(labelsize=0)
# ax.set_ylabel('Height [km]')
# ax.set_title(u'{0} {1}\N{DEGREE SIGN} azimuth cross-section \n{2}'.
#              format(pyart.graph.common.generate_radar_name(radar), azimuth,
#                     pyart.util.datetime_from_radar(xsect).strftime('%d %b %Y %H:%M:%S UTC')), loc='left')
# ax.set_title('Reflectivity [dBZ]', loc='right')
#
# plt.show()
