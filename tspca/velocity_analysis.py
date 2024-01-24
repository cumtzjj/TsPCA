import obspy
import numpy as np
from obspy.core import Trace
from scipy.signal import hilbert
from scipy.interpolate import interp1d


"""
velocity analysis for autocorrelation data
"""
def velocity_analysis(ac_data, ray_para, vp_num=401, vp_max=10, PWS_order=0):
    # calculate the number of autocorrelation data
    ac_event_num = len(ac_data)        
    # calculate the length of single autocorrelation data
    ac_single_num = ac_data[0].stats.npts
    # variables for velocity analysis
    vp = np.linspace(0, vp_max, vp_num)
    # velocity spectrum
    ac_nmo = np.zeros([ac_single_num, ac_event_num])
    ac_vel_spec = np.zeros([ac_single_num, vp_num])
    for i in range(vp_num):
        vp_temp = vp[i]
        for j in range(ac_event_num):
            t = ac_data[j].times()
            rayp_temp = ray_para[j]
            # calculate arrival time for each ray parameter
            tt = t * np.sqrt(1 - rayp_temp**2 * vp_temp**2)
            # interpolate data
            ac_nmo[:, j] = interp1d(t, ac_data[j].data, kind='cubic')(tt)
            # calculate average amplitude
        # stack velocity spectrum
        if PWS_order == 0:
            ac_vel_spec[:, i] = np.mean(ac_nmo, axis=1)
        else:
            stack = np.zeros(ac_single_num)
            phase = 0j
            for j in range(ac_event_num):
                stack += ac_nmo[:, j]
                # calculate phase
                asig = hilbert(ac_nmo[:, j])
                phase += asig / np.abs(asig)
            stack /= ac_event_num
            weight = np.abs((phase) / ac_event_num)
            ac_vel_spec[:, i] = stack*weight**PWS_order
    return ac_vel_spec


"""
moveout correction for autocorrelation data
"""
def moveout_corr(ac_data, ray_para, Va):
    # calculate the number of autocorrelation data
    ac_event_num = len(ac_data)        
    # corrected autocorrelation stream
    ac_stream = obspy.Stream()
    for j in range(ac_event_num):
        t = ac_data[j].times()
        rayp_temp = ray_para[j]
        # calculate corrected arrival time
        tt = t * np.sqrt(1 - rayp_temp**2 * Va**2)
        # interpolate data
        ac_nmo = interp1d(t, ac_data[j].data, kind='cubic')(tt)
        ac_stream.append(Trace(header=ac_data[j].stats, data=ac_nmo))
    return ac_stream