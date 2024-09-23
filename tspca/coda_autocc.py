import numpy as np
from obspy.signal.util import next_pow_2
from scipy import signal
from obspy.core import Stream


def comp_autocc(tr,W = 0.1,gauss = 3, freq_min = 0.1, freq_max = 2.0, out_length = 40):
    """
    Computing the autocorrelation function with spectral whitenning
    """

    centre_gauss = 0
    t_delay = 10

    fs = 1/(tr.stats.sac.delta)
    fs = np.round(fs)
#    print(fs)
    fs2 = 0.5*fs
    b, a = signal.butter(4, [freq_min/fs2,freq_max/fs2], 'bandpass')

    data = tr.data
    
    m = tr.stats.npts
#    print(m)
    wt = signal.windows.tukey(m,0.2)
    data = data*wt

    n = next_pow_2(m)
    n *= 2

#    print(n)
    data = np.fft.fft(data,n)
    conz = np.conj(data)
    data = data*conz

    df = 1/(n*tr.stats.sac.delta)
    dw = 2*np.pi*df
    w = np.arange(n)*dw
    a2 = 4*gauss**2


    gw = np.exp(-((w-centre_gauss*np.pi*2)**2)/a2-w*t_delay*1j)
    data2 = data*gw

    N=np.int(W/(2*df))
#    print(N)

    d=np.zeros(n)
    for i in range(n):
        N0=np.min([n,i+N])-np.max([0,i-N])
        d[i]=1/N0*np.sum(np.abs(data[np.max([0,i-N]):np.min([n,i+N])]))

#    if any(d < 1e-8): raise Exception('Zero division')

    data_dec = data2/d
    data_dec[d<1e-8] = 0

    data_time = np.fft.ifft(data_dec).real

    data_time = data_time[:m]/np.max(data_time)

    data_taper = data_time*wt
    
    data_bp = signal.filtfilt(b,a,data_taper)


    tr.data = -data_bp[np.int(t_delay*fs):np.int(t_delay*fs+out_length*fs+1)]


 #   print(len(tr.data))

    return tr
