# TsPCA
![LICENSE](https://img.shields.io/badge/license-MIT-green)

data processing for Teleseismic P-wave Coda Autocorrelation
made by Tianyu Cui and Jinju Zhou, January 24, 2024

If you intent to use a part or entire code in your research, please cite at least one of the following articles:
Zhou, J., & Zhang, W. (2021). Extracting reliable P-wave reflections from teleseismic P wave coda autocorrelation. Journal of Geophysical Research: Solid Earth, 126, e2021JB022064. https://doi.org/10.1029/2021JB022064
Zhou, J., Hu, N., Hu, Y., & Zhang, W. (2022). Feasibility of estimating parameters and enhancing reflections of dipping Moho from teleseismic P wave coda autocorrelation. Journal of Geophysical Research: Solid Earth, 127, e2021JB023688. https://doi.org/10.1029/2021JB023688

## features:
1. spectral whitening and autocorrelation for teleseismic P-wave coda wave 
2. velocity analysis to find average velocity above interface and two-way traveltimes
3. moveout corrected for teleseismic P-wave coda autocorrelation


## dependencies
obspy, signal, numpy, pandas, matplotlib, pygmt
***

## parameters setting
```python

    ### parameters setting
    data_wave_folder = "demo_data"   # folder of waveforms
    # data preprocessing setting
    channel = 'Z'           # channel of waveforms
    cut_win = [-20, 80]     # cut window for autocorrelation, before 20s and after 80s of P wave. unit: second
    signal_win = [-5, 5]    # signal window for calculating SNR, 5s before and after P wave. unit: second
    noise_win = [-20, -5]   # noise window for calculating SNR, 20s before and 5s before P wave. unit: second
    freq_pre = [0.1, 2]     # frequency band for prefiltering data, 0.1-2Hz
    min_SNR = 10            # minimum SNR for autocorrelation
    # spectral whitening and autocorrelation setting
    dw = 0.3                # spectral whitening width
    filter_data = True      # filter autocorrelation data or not
    freq_ac = [0.2, 2]      # frequency band for autocorrelation data, 0.1-1Hz
    taper_data = True       # taper autocorrelation data or not
    taper_length = 5        # taper length for autocorrelation data, 5s
    data_PWS_order = 1      # order of phase-weighted stack for autocorrelation data
    # velocity analysis setting
    vel_PWS_order = 1       # order of phase-weighted stack for velocity analysis data, 0 for linear stack
    vp_max = 10             # maximum velocity for velocity analysis, default: 10km/s
    vp_scan_num = 401               # number of velocity for velocity analysis
    max_scan_time_win = [10, 15]    # time window for scanning maximum velocity and t0, 10-15s of autocorrelation data
    # plot figure setting
    figure_folder = 'figure'    # data folder for saving figures
    ac_view_cut = 20            # time window for plotting autocorrelation and velocity spectrum, 20s of autocorrelation data
    plot_sta_event = True       # plot station and event distribution or not

```

## an example
```bash
    python coda_process.py
```
