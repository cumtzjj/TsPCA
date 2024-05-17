import os
import glob
import obspy
import numpy as np
from obspy.taup import TauPyModel

"""
convert earthquake data to station centered data
cut data from referenced arrival time and filter data 
"""
def event2station(data_dir):
    # list all same station file
    data_file = [data for data in glob.glob('%s/*/*.sac'%(data_dir))]
    # list all station names
    single_sta = sorted(list(set([obspy.read(data_sac)[0].stats.station for data_sac in data_file])))
    # loop for each station
    for sta in single_sta:
        # mkdir for each station
        if not os.path.exists('wave_station'):
            os.makedirs('wave_station')
        if not os.path.exists('wave_station/%s'%(sta)):
            os.makedirs('wave_station/%s'%(sta))
        # list all sac file for each station
        sacfiles_sta = sorted(glob.glob('%s/*/*.%s.*sac'%(data_dir,sta)))
        # copy sac file
        for sac in sacfiles_sta:
            os.system('cp %s %s'%(sac, 'wave_station/%s/'%(sta)))


"""
calculate signal to noise ratio (SNR)
"""
def SNR_calculate(data, ref_time, signal_win=[-5, 5], noise_win=[-20, 0]):
    # calculate SNR
    signal = data.copy()
    noise = data.copy()
    # cut signal and noise
    starttime = data.stats.starttime
    signal.trim(starttime+ref_time+signal_win[0], starttime+ref_time+signal_win[1])
    noise.trim(starttime+ref_time+noise_win[0], starttime+ref_time+noise_win[1])
    # calculate SNR
    noise_power = np.sum(noise.data**2)
    signal_power = np.sum(signal.data**2)
    snr = 10 * np.log10(signal_power / noise_power)
    return snr


'''
return predicted phase arrival time 
'''
def predicted_phase_arrival(event_depth, distance, phase_list, ref_model="ak135"):
    model = TauPyModel(model=ref_model)
    arrival_phase = {}
    arrivals = model.get_travel_times(source_depth_in_km=event_depth, distance_in_degree=distance, phase_list=phase_list)
    for index, arrival in enumerate(arrivals):
        phase_name = arrival.name
        # choose the first arrival time of P arrival
        if phase_name == "P" and index == 0:
            phase_arrival = arrival.time
            arrival_phase[phase_name] = phase_arrival
        elif phase_name == "P" and index != 0:
            continue
        else:
            phase_arrival = arrival.time
            arrival_phase[phase_name] = phase_arrival
    return arrival_phase


"""
cut data from reference P arrival time and filter data
select data with signal-to-ratio (SNR)
"""
def data_select(data_trace, cut_win=[-20, 60], signal_win=[-5, 5], noise_win=[-20,-5], filter_min=0.1, filter_max=2, filter_corner=4):
    # copy data
    tr = data_trace.copy()
    # detrend data
    tr.detrend('demean')
    tr.detrend('linear')
    # resample data
    tr.resample(40)
    # filter data
    tr.filter('bandpass', freqmin=filter_min, freqmax=filter_max, corners=filter_corner, zerophase=True)
    # event depth and distance (degree)
    evdp = tr.stats.sac.evdp
    dist = tr.stats.sac.gcarc
    # predicted P arrival time
    ref_P_time = predicted_phase_arrival(evdp, dist, ['P'])['P']
    # cut data from reference P arrival time
    tr.trim(tr.stats.starttime+ref_P_time+cut_win[0], tr.stats.starttime+ref_P_time+cut_win[1], pad=True, fill_value=0)
    # calculate SNR
    snr = SNR_calculate(data_trace.copy(), ref_P_time, signal_win, noise_win)
    return tr, snr


"""
return maximum velocity and two-way arrival time in time window
"""
def get_max_value(vel_spec_data, max_scan_time_win, vp_max, time_max):
    print('   scanning maximum velocity and t0 in time window %.1f-%.1f s'%(max_scan_time_win[0], max_scan_time_win[1]))
    min_scan_time = max_scan_time_win[0]
    max_scan_time = max_scan_time_win[1]
    vp_ori_num = vel_spec_data.shape[0]
    t0_ori_num = vel_spec_data.shape[1]
    min_scan_time_index = int(min_scan_time / time_max * t0_ori_num)
    max_scan_time_index = int(max_scan_time / time_max * t0_ori_num)
    vel_data = vel_spec_data[min_scan_time_index:max_scan_time_index, :]
    vel_data = vel_spec_data[:, min_scan_time_index:max_scan_time_index]
    max_Va_index, max_t0_index = np.unravel_index(np.argmax(vel_data), vel_data.shape)
    vp_num, t0_num = vel_data.shape
    max_Va = np.linspace(0, vp_max, vp_num)[max_Va_index]
    max_t0 = np.linspace(min_scan_time, max_scan_time, t0_num)[max_t0_index]
    print('   maximum average velocity (Va): %.2f km/s, two-way time (t0): %.2f s'%(max_Va, max_t0))
    return max_Va, max_t0