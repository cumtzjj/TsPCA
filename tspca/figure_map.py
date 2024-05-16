import os
import numpy as np
import matplotlib.pyplot as plt

"""
plot station and events distribution
"""
def station_center_map(stname, min_SNR, stla, stlo, evla, evlo, figure_folder='figure'):
    try:
        import pygmt
    except:
        raise ImportError("pygmt is not installed, please install it first.")
    print('   plot station and events distribution map...')
    figure_path = os.path.join(figure_folder)
    if not os.path.exists(figure_folder):
        os.makedirs(figure_path)
    pygmt.config(FONT_TITLE="10p,Helvetica,black")
    fig = pygmt.Figure()
    fig.coast(projection="E%s/%s/95/10c"%(stlo, stla), region="g", land="grey", water="white", area_thresh=10000, 
              frame=["f", "+t%s: total %s events (min SNR: %s)"%(stname, len(evla),min_SNR)])
    # map distance range 30, 60 and 90 degree
    for style in ["E-60d", "E-120d", "E-180d"]:
        fig.plot(x=stlo, y=stla, style=style, pen="0.8p,red" )
    fig.plot(x=evlo, y=evla, style="a0.25c", fill="red", pen="black")
    fig.plot(x=stlo, y=stla, style="t0.3c", fill="black", pen="black")
    # text 30, 60 and 90 degree
    for y,text in zip([14, -16, -46], ["30", "60", "90"]):
        fig.text(x=stlo, y=y, text=text, offset="0c/0.1c")
    fig.savefig('%s/station_events_%s.png'%(figure_path,stname))

"""
plot velocity spectrum figure
"""
def plot_vel_spec(ac_data, vel_spec_data, max_Va, max_t0, vp_max=10, time_cut_max=20, figure_folder='figure', ):
    # figure for each station
    figure_path = os.path.join(figure_folder)
    if not os.path.exists(figure_folder):
        os.makedirs(figure_path)
    sta_name = ac_data[0].stats.station
    ac_max_time = np.max(ac_data[0].times())
    # plot velocity spectrum
    plt.figure(figsize=(6, 4))
    vel_data = vel_spec_data.T
    vmin, vmax = -np.max(np.abs(vel_spec_data)), np.max(np.abs(vel_spec_data))
    im = plt.imshow(vel_spec_data.T, aspect='auto', cmap='jet', extent=[0, ac_max_time, vp_max, 0], vmin=vmin, vmax=vmax)
    plt.text(max_t0, max_Va, '%.2f, %.2f' % (max_t0, max_Va), color='black', fontsize=10)
    plt.scatter(max_t0, max_Va, s=30, marker='+', color='cyan', lw=1.2,)
    plt.xlim(0, time_cut_max)
    plt.ylim(vp_max, 0)
    plt.xlabel('t$_0$ (s)')
    plt.ylabel('V$_a$ (km/s)')
    plt.colorbar(im)
    plt.savefig('%s/vel_%s.png' % (figure_path, sta_name))
    plt.close()

"""
plot autocorrelation and PWS figure
"""
def plot_ac(ac_data, stack_data, time_cut_max=20, taper_length=5, taper=False, figure_folder='figure', moveout_corr=False, ):    
    # figure for each station
    figure_path = os.path.join(figure_folder)
    if not os.path.exists(figure_folder):
        os.makedirs(figure_path)
    sta_name = ac_data[0].stats.station        
    # plot autocorrelation and PWS
    stack = stack_data
    time = np.arange(stack.stats.npts) * stack.stats.delta
    fig, (ax_auto, ax_pws) = plt.subplots(1, 2, figsize=(6, 4), gridspec_kw={'width_ratios': [4, 1]})
    ax_auto.set_xlim(29, 95)
    ax_auto.set_ylim(time_cut_max, 0)
    ax_auto.set_ylabel('Time (s)')
    ax_auto.set_xlabel('Epicentral distance ($^\circ$)')
    ax_auto.set_title('Vertical autocorrelograms (%s station)'%sta_name)
    for auto in ac_data:
        auto.data /= np.max(np.abs(auto.data))
        epi_dist = auto.stats.sac.gcarc
        ax_auto.plot(auto.data + epi_dist, time, lw=0.5, color='black')
        ax_auto.fill_betweenx(time, epi_dist, auto.data + epi_dist, lw=0.5, color='gray', where=(auto.data > 0))
    ax_pws.set_ylim(time_cut_max, 0)
    ax_pws.set_xlim(-1, 1)
    ax_pws.set_title('PWS')
    stack.data /= np.max(np.abs(stack.data))
    ax_pws.plot(stack.data, time, lw=0.5, color='black')
    ax_pws.fill_betweenx(time, 0, stack.data, lw=0.5, color='gray', where=(stack.data > 0))
    if moveout_corr:
        plt.savefig('%s/ac_%s_moveout.png'%(figure_path, sta_name))
    else:
        plt.savefig('%s/ac_%s.png'%(figure_path, sta_name))
    plt.close()