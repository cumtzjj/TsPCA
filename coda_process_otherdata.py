#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import obspy
import numpy as np
from obspy.taup import TauPyModel
import tspca

"""
calculate the autocorrelation of each station and plot cross-correlation figures
"""
def coda_autocorrelation(datafolder, channel='Z', cut_win=[-20, 80], signal_win=[-5, 5], noise_win=[-20, -5], min_SNR=3, ac_view_cut=20, taper_length=5,
                         freq_pre=[0.1, 2], dw=0.5, data_PWS_order=1, vel_PWS_order=0, freq_ac=[0.2, 2], vp_max=10, vp_scan_num=401, filter_data=True, 
                         taper_data=False, max_scan_time_win=[10, 15], plot_sta_event=True, out_length = 40, gauss = 5, figure_folder='figure'):    
    # MAIN PROCESSING FOR VERTICAL AUTOCORRELOGRAMS
    valid_station = sorted(list(set([obspy.read(sacdata)[0].stats.station for sacdata in glob.glob('%s/*/*.SAC'%(datafolder))])))
    print(valid_station)
    # information for each station
    auto_all_stations = []
    auto_corr_all_stations = []
    stack_all_stations = obspy.Stream()
    stack_corr_all_stations = obspy.Stream()
    vel_spec_all_stations = []
    station_elevation = []
    max_Va_list = []
    max_t0_list = []
    model = TauPyModel(model="iasp91")
    for station in valid_station: 
        print(station)
#        data_file = os.path.join(datafolder, station)
        # autocorrelation stream 
        auto_stream = obspy.Stream()
        # autocorrelation stream for velocity analysis
        auto_stream_vel = obspy.Stream()
        # ray parameter
        ray_parameter = []
        # selected data
        selected_data = []
        # station and event coordinates
        sta_event_coord = {'evla': [], 'evlo': []}
        # select channel data
#        sta_data_files = glob.glob('%s/*.sac'%(data_file))
        sta_data_files = glob.glob('%s/*/*%s*.SAC'%(datafolder,station))
        sta_data_files = [data for data in sta_data_files if obspy.read(data)[0].stats.channel[-1] == channel]
        for fname in sorted(sta_data_files):
            print(fname)
            tr = obspy.read(fname, format='SAC')[0]
            print('delta = %f'%tr.stats.delta)
            #print('npts = %d'%tr.stats.npts)
            # delete data with no data or nan or inf
            if np.any(tr.data) == 0 or np.isnan(tr).any() or np.isinf(tr).any():
                print('%s has no data and delet!'%fname)
                os.system('rm %s'%fname)
            if hasattr(tr.stats.sac,'evdp'):
                print(tr.stats.sac.evdp)
            # delete data with low SNR
                tr_cut, SNR, SNR2 = tspca.data_select(tr, cut_win, signal_win, noise_win, freq_pre[0], freq_pre[1], 4)
                print(SNR, SNR2)
                if SNR > min_SNR and SNR2 > 2:
                # station and event coordinates
                    sta_event_coord['stla'] = tr.stats.sac.stla
                    sta_event_coord['stlo'] = tr.stats.sac.stlo
                    sta_event_coord['stel'] = tr.stats.sac.stel
                    sta_event_coord['evla'].append(tr.stats.sac.evla)
                    sta_event_coord['evlo'].append(tr.stats.sac.evlo)
                # calculate ray parameter
                    arrival_info = model.get_travel_times(source_depth_in_km=tr.stats.sac.evdp, distance_in_degree=tr.stats.sac.gcarc, phase_list=['P'])[0]
                    ray_para = arrival_info.ray_param / 6371
                    ray_parameter.append(ray_para)
                # compute velocity analysis
                #auto_vel = tspca.compute_auto(tr_cut.copy(), dw, freq_ac[0], freq_ac[1], taper_length, taper_data, filter_data, vel_analysis=True)
                #auto_stream_vel.append(auto_vel)
                # compute auto-correlation
                    auto = tspca.comp_autocc(tr_cut.copy(), dw, gauss, freq_ac[0], freq_ac[1], out_length)
                    print(auto)
                    auto_stream.append(auto)
                    selected_data.append(tr_cut.copy().data)
            else: print('evdp not exist')
        np.save('selected_data_%s.npy'%station,selected_data)
        # station and events plot
        print('1. auto-selected events (SNR>= %.1f): %d/%d'%(min_SNR, len(auto_stream), len(sta_data_files)))  
        print('2. calculate autocorrelation for %s station'%station)
        if plot_sta_event:
            tspca.station_center_map(station, min_SNR, sta_event_coord['stla'], sta_event_coord['stlo'], sta_event_coord['evla'], sta_event_coord['evlo'])
        # velocity analysis
        print('3. velocity analysis for %s station'%station)
        vel_spec = tspca.velocity_analysis(auto_stream, ray_parameter, vp_num=vp_scan_num, vp_max=vp_max, PWS_order=vel_PWS_order)
        ac_time_len = out_length
        # get max Va and t0
        max_Va, max_t0 = tspca.get_max_value(vel_spec.T, max_scan_time_win, vp_max, ac_time_len)
        max_Va_list.append(max_Va)
        max_t0_list.append(max_t0)
        # moveout correction
        print('4. moveout correction for %s station, based on Va=%.2f km/s'%(station, max_Va))
        auto_corr = tspca.moveout_corr(auto_stream, ray_parameter, max_Va)
        stack_corr = tspca.data_stack(auto_corr, data_PWS_order)
        auto_corr_all_stations.append(auto_corr)
        stack_corr_all_stations.append(stack_corr)
        vel_spec_all_stations.append(vel_spec)
        auto_all_stations.append(auto_stream)
        stack = tspca.data_stack(auto_stream, data_PWS_order)
        stack_all_stations.append(stack)
        # station elevation
        station_elevation.append(auto.stats.sac.stel/1000)
        # plot autocorrelation and velocity spectrum
        print('5. plot autocorrelation, velocity spectrum and moveout corrected autocorrelation...')
        tspca.plot_vel_spec(auto_stream, vel_spec, max_Va, max_t0, vp_max, ac_view_cut, figure_folder)
        tspca.plot_ac(auto_stream, stack, ac_view_cut, taper_length, taper_data, figure_folder)
        tspca.plot_ac(auto_corr, stack_corr, ac_view_cut, taper_length, taper_data, figure_folder, moveout_corr=True)
        # binned stack
        binned_stack(auto_stream,8,station)
    # save autocorrelation and velocity spectrum
    print('6. save autocorrelation and velocity spectrum...')
    data_total_info = {'stations': valid_station, 'auto_all_stations': auto_all_stations, 'auto_corr_all_stations': auto_corr_all_stations, 'stack_all_stations': stack_all_stations,
                       'stack_corr_all_stations': stack_corr_all_stations, 'vel_spec_all_stations': vel_spec_all_stations, 'station_elevation': station_elevation,
                       'max_Va_list': max_Va_list, 'max_t0_list': max_t0_list, }
    np.save('data_total_info.npy', data_total_info)


def binned_stack(tr, nbins,station):
    """

    """
#    ev_depth=tr.stats.sac.evdp
    npts=tr[0].stats.npts
    model = TauPyModel(model="iasp91")
    rayp=np.zeros([len(tr),1])
    obs=np.zeros([npts,len(tr)])
    baz=np.zeros([len(tr),1])
    k=0
    for auto in tr:
        ev_depth=auto.stats.sac.evdp
        epi_dist=auto.stats.sac.gcarc

        arrivals = model.get_travel_times(source_depth_in_km=ev_depth,
            distance_in_degree=epi_dist,
            phase_list=["P"])

        arrs=arrivals[0]


        rad2km=6371
        rayp[k]=arrs.ray_param/rad2km
        obs[:,k]=auto.data
        baz[k]=auto.stats.sac.baz
        k+=1


    bins = np.linspace(0.04,0.08,nbins+1)

    print(bins)

    inds = np.digitize(rayp,bins)

    print(inds)

    zac_bins=np.zeros([npts,nbins])
    zac_std=np.zeros([npts,nbins])
    rayp_bins=np.zeros(nbins)
    events_num=np.zeros(nbins)

    for j in range(nbins):
        inx=[i for i,x in enumerate(inds) if x==j+1]
        print(len(inx))
        if len(inx)>0:
            rayp_bins[j]=np.mean(rayp[inx])
            events_num[j]=len(inx)
            ac_bins=obs[:,inx]
            baz_bins=baz[inx]
            slow_bins=rayp[inx]
            np.save('baz_bins_%s_%s'%(str(j),station),arr=baz_bins)
            np.save('zac_bins_%s_%s'%(str(j),station),arr=ac_bins)
            np.save('slow_bins_%s_%s'%(str(j),station),arr=slow_bins)
            zac_bins[:,j]=np.mean(ac_bins,axis=1)
            zac_std[:,j]=np.std(ac_bins,axis=1)

    np.save('zac_bins_%s.npy'%station,arr=zac_bins)
    np.save('zac_std_%s.npy'%station,arr=zac_std)
    np.save('rayp_bins_%s.npy'%station,arr=rayp_bins)
    np.save('events_num_%s.npy'%station,arr=events_num)

if __name__ == '__main__':

    ### parameters setting
    data_wave_folder = "/mnt/d/work/DownloadSeisData-main/Bash/DATA/SAC_event/II"   # folder of waveforms
    # data preprocessing setting
    channel = 'Z'           # channel of waveforms
    cut_win = [-20, 80]     # cut window for autocorrelation, before 20s and after 80s of P wave. unit: second
    signal_win = [-5, 5]    # signal window for calculating SNR, 5s before and after P wave. unit: second
    noise_win = [-20, -10]   # noise window for calculating SNR, 20s before and 5s before P wave. unit: second
    freq_pre = [0.1, 2]     # frequency band for prefiltering data, 0.1-2Hz
    min_SNR = 10             # minimum SNR for autocorrelation
    # spectral whitening and autocorrelation setting
    dw = 0.1                # spectral whitening width
    filter_data = True      # filter autocorrelation data or not
    freq_ac = [0.2, 1]      # frequency band for autocorrelation data, 0.1-1Hz
    taper_data = True       # taper autocorrelation data or not
    taper_length = 5        # taper length for autocorrelation data, 5s
    data_PWS_order = 1      # order of phase-weighted stack for autocorrelation data
#    out_length = 40         # length of output autocorrelation 
    # velocity analysis setting
    vel_PWS_order = 1       # order of phase-weighted stack for velocity analysis data, 0 for linear stack
    vp_max = 10             # maximum velocity for velocity analysis, default: 10km/s
    vp_scan_num = 401               # number of velocity for velocity analysis
    max_scan_time_win = [10, 20]    # time window for scanning maximum velocity and t0, 10-15s of autocorrelation data
    # plot figure setting
    figure_folder = 'figure'    # data folder for saving figures
    ac_view_cut = 20            # time window for plotting autocorrelation and velocity spectrum, 20s of autocorrelation data
    plot_sta_event = False #True       # plot station and event distribution or not
    out_length = 40
    gauss = 3


    ### output:
    # 1. autocorrelation and velocity spectrum figures for each station (under figure folder)
    # 2. autocorrelation and velocity spectrum data for each station (data_total_info.npy)

    # calculate autocorrelation and velocity analysis
    coda_autocorrelation(data_wave_folder, channel, cut_win, signal_win, noise_win, min_SNR, ac_view_cut, taper_length,
                         freq_pre, dw, data_PWS_order, vel_PWS_order, freq_ac, vp_max, vp_scan_num, filter_data, 
                         taper_data, max_scan_time_win, plot_sta_event, out_length, gauss, figure_folder)
    
