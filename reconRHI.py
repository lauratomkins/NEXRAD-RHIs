#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 11:11:55 2021

@author: ltomkins
"""

from utils import rhi_fun
import pyart
import matplotlib.pyplot as plt
import cmocean
import matplotlib.cm as cm

#%%
# radar object (PPI)
radar = pyart.io.read_nexrad_archive("data/KBGM20200207_132642_V06")

# RHI object
xsect = rhi_fun.constructRHI(radar, 120, 'reflectivity')

# Plot
display_ppi = pyart.graph.RadarDisplay(radar)
display_rhi = pyart.graph.RadarDisplay(xsect)

#%%
fig = plt.figure(figsize=(10, 3))
ax = plt.axes()

display_rhi.plot('reflectivity', vmin=5, vmax=45, colorbar_flag=False, title_flag=False,
                 cmap=cmocean.tools.crop_by_percent(cm.magma_r, 10, which='max', N=None))
display_rhi.set_aspect_ratio(5)
display_rhi.set_limits(xlim=(0, 200), ylim=(0, 6))
ax.set_ylabel('Height [km]')
ax.set_title(u'{0} {1}\N{DEGREE SIGN} Azimuth Cross Section \n{2}'.
             format(pyart.graph.common.generate_radar_name(radar), str(120),
                    pyart.util.datetime_from_radar(xsect).strftime('%d %b %Y %H:%M:%S UTC')), loc='left')
ax.set_title('Reflectivity [dBZ]', loc='right')

plt.show()
