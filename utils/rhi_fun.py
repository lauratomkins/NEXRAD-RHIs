import pyart
from utils import rad_fun

def constructRHI(radar, azimuth, field):

    # Find tilts with data for given field
    all_tilts, all_angs = rad_fun.returnAll(radar, field)

    # Due to the scan strategy, we have to remove the duplicate tilts
    good_tilts = rad_fun.removeDuplicate(all_angs, all_tilts)
    radar = radar.extract_sweeps(good_tilts)

    # generate cross section (RHI)
    xsect = pyart.util.cross_section_ppi(radar, [azimuth])

    return xsect