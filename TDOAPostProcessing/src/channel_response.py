import numpy as np
from scipy.signal import group_delay
from scipy import fftpack
import scipy
from src.wave2toa import recwav
from matplotlib import pyplot as plt

if __name__ == '__main__':
    matlab_path = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/allchirps_channel_measurement.mat'
    record = recwav()
    record.change_path('/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/channel_estimation.WAV', 'in')
    record.PlotSignal()


    #load matlab expected signal
    matlab_chirps = scipy.io.loadmat(matlab_path)
    chirps = matlab_chirps['allchirp']

    spectf, spectt, spectsxx = record.Spectogram()
    # cut it by the designed record:
    channels = {}
    signals = []
    sec_in_taps = record.sample_rate
    for i in [1,2,3,4]:
        signals.append(record.signal[(i-1) * sec_in_taps: i * sec_in_taps + 1])
        length_of_fft = len(chirps[0]) + len(signals[0]) + 1
        mat_fft = fftpack.fft(chirps[i-1], length_of_fft)
        sig_fft = fftpack.fft(signals[i-1], length_of_fft)
        channel_fft = sig_fft / mat_fft
        channel_time = scipy.ifft(channel_fft)
        # calculate Group delay
        gd = group_delay((channel_time, 1), 60000)
        channels['sp' + str(i)] = {'matlab_time': chirps[i-1],
                                   'matlab_frequency': mat_fft,
                                   'record_time': signals[i-1],
                                   'record_frequency': sig_fft,
                                   'channel_time': channel_time,
                                   'channel_frequency': fftpack.fftshift(channel_fft),
                                   'group_delay': gd
                                   }



        fig = plt.figure(1)
        fig2 = plt.figure(2)
        ax = fig.add_subplot(1,1,1)
        ax.plot(channels['sp1']['channel_frequency'])
        ax.plot(channels['sp2']['channel_frequency'])
        ax.plot(channels['sp3']['channel_frequency'])
        ax.plot(channels['sp4']['channel_frequency'])
        ax.title('Channels Group Delay')
        ax.legend(['sp1', 'sp2', 'sp3', 'sp4'])
        ax.grid
        ax2 = fig2.add_subplot(1,1,1)
        ax2.plot(channels['sp1']['group_delay'])
        ax2.plot(channels['sp2']['group_delay'])
        ax2.plot(channels['sp3']['group_delay'])
        ax2.plot(channels['sp4']['group_delay'])
        ax2.legend(['sp1', 'sp2', 'sp3', 'sp4'])
        ax2.title('Channels Group Delay')
        ax2.grid

        fig2.show()
        fig.show()


    pass





