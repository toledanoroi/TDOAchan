import csv
import wave
import numpy as np
from pandas import DataFrame as df
import time
import os
from scipy import signal
from scipy import io
import scipy.fftpack as fftpkt
from matplotlib import pyplot as plt
import sys
from termcolor import colored
import src.algorithms as algos
from utils import UTILS
from utils import SignalHandler
from collections import OrderedDict as OD

plotting = False

def resampling(speaker,rec_Fs):
    factor = int((rec_Fs/250000.0)*len(speaker.matlab_chirp))
    f = signal.resample(speaker.matlab_chirp,factor)
    x = np.linspace(0, 0.001, 251, endpoint=False)
    x1 = np.linspace(0, 0.001, factor, endpoint=False)
    if plotting:
        plt.plot(x, speaker.matlab_chirp, 'go-', x1, f, '.-',x1,speaker.chirp,'.-', 1, speaker.matlab_chirp[0], 'ro')
        plt.legend(['data', 'resampled','python generated'], loc='best')
        plt.show()

    fft1   = fftpkt.fftshift(fftpkt.fft(speaker.matlab_chirp))
    freqs  = fftpkt.fftshift(fftpkt.fftfreq(len(speaker.matlab_chirp), 1.0 / rec_Fs))
    fft22  = fftpkt.fftshift(fftpkt.fft(f))
    freqs2 = fftpkt.fftshift(fftpkt.fftfreq(len(f), 1.0 / rec_Fs))
    fft33  = fftpkt.fftshift(fftpkt.fft(speaker.chirp))
    freqs3 = fftpkt.fftshift(fftpkt.fftfreq(len(speaker.chirp), 1.0 / rec_Fs))

    if plotting:
        plt.figure(10)
        plt.plot(freqs, abs(fft1), 'r')
        plt.plot(freqs2, abs(fft22), 'b')
        plt.plot(freqs3, abs(fft33), 'g')
        plt.legend(['matlab 250000', 'resampled', 'python generated'], loc='best')
        plt.show()


class recwav(object):
    def __init__(self):
        self.path = '../inputs/blackmanharris5ms/1.wav'
        self.sample_rate = 96000
        self.toa_csv_path = '../output/toa_record_'+str(int(time.time()))+'.csv'
        self.results_path = '../output/locations_results_'+str(int(time.time()))+'.csv'
        self.record_sig = []
        self.num_of_channels = 1
        self.prev_time = 0
        self.curr_time = 0
        self.parttime = 0.49
        self.total_time = 0
        self.time_samples_vect = []

    def PlotSignal(self,title):
        self.sig_file_pointer = wave.open(self.path,'r')
        self.num_of_channels = self.sig_file_pointer.getnchannels()
        # Extract Raw Audio from Wav File
        self.signal = np.fromstring(self.sig_file_pointer.readframes(-1), 'Int16')
        self.sample_rate = self.sig_file_pointer.getframerate()

        # If Stereo
        if self.num_of_channels == 2:
            print 'Our Application support only mono files , record from one input only'
            sys.exit(0)

        self.rec_time = np.linspace(0, len(self.signal) / self.sample_rate , num=len(self.signal))
        self.total_time = self.rec_time[-1]

        print "Signal Parameters:\n\tFile Name:\t{0}\n\tSample rate:\t{1}" \
              "\n\tNumber of channels:\t{2}\n\tDuration Time:\t{3}".format(self.path,
                                                                           self.sample_rate,
                                                                           self.num_of_channels,
                                                                           self.total_time)

        if plotting:
            plt.figure(1)
            plt.title('Signal Wave: ' + title )
            plt.plot(self.rec_time, self.signal)
            plt.show()

    def PlotPartSig(self,title):

        # If Stereo
        if self.num_of_channels == 2:
            print 'Our Application support only mono files'
            sys.exit(0)

        self.part_rec_time = np.linspace(0, float(len(self.curr_sig)) / self.sample_rate, num=len(self.curr_sig))

        plt.figure(12)
        plt.title('Signal Wave: ' + title )
        plt.plot(self.part_rec_time,self.curr_sig)
        plt.show()

    def PlotFFT(self,title):
        self.fft = fftpkt.fft(self.signal)
        self.psd = np.abs(self.fft) ** 2
        self.freqs = fftpkt.fftshift(fftpkt.fftfreq(len(self.signal),1.0/self.sample_rate))
        self.fft = fftpkt.fftshift(self.fft)

        if plotting:
            plt.figure()
            plt.title('Signal Wave FFT: ' + title)
            plt.plot(self.freqs, abs(self.fft))
            plt.grid()
            plt.show()

    def change_path(self, path, mode='in'):
        if mode == 'in':
            self.path = path
        else:
            self.toa_csv_path = path

    def plotOnlyOneCycle(self):
        tmpsig = self.signal[int(0.5*self.sample_rate):int(1.3*self.sample_rate)]
        tmpfft = fftpkt.fftshift(fftpkt.fft(tmpsig))
        tmppsd = np.abs(tmpfft) ** 2
        tmpfreqs = fftpkt.fftshift(fftpkt.fftfreq(len(tmpsig), 1.0 / self.sample_rate))
        tmprec_time = np.linspace(0, len(tmpsig) / self.sample_rate, num=len(tmpsig))

        if plotting:
            plt.figure(3)
            plt.title('one cycle')
            plt.plot(tmprec_time,tmpsig)
            plt.show()

            plt.figure(4)
            plt.title('one cycle FFT:')
            plt.plot(tmpfreqs,abs(tmpfft))
            plt.grid()
            plt.show()

    def Spectogram(self):
        f,t,Sxx = signal.spectrogram(self.signal,self.sample_rate, nperseg=64)
        plt.pcolormesh(t,f,(10 * np.log10(Sxx)), cmap='inferno')
        plt.colorbar()
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        plt.show()

    def CutCurrSig(self):
        self.total_time = self.rec_time[-1]
        self.curr_time += self.parttime
        self.time_samples_vect.append(self.curr_time)
        if self.curr_time > self.total_time:
            print "curr_time exceed total time , update part to total time"
            self.curr_time = self.total_time
        self.curr_sig = self.signal[int(self.prev_time*self.sample_rate): int(self.curr_time*self.sample_rate)]
        self.PlotPartSig(str(self.curr_time))
        self.prev_time = self.curr_time
        return self.curr_sig

    def CutSigByPeaks(self, ind, peaks, filter_size):
        if ind >= (len(peaks) - 1):
            self.curr_sig = self.signal[peaks[ind] - filter_size:]
        else:
            final = peaks[ind + 1] - 500
            self.curr_sig = self.signal[peaks[ind] - filter_size : final]
        # calculate FFT
        self.curr_sig_fft = fftpkt.fftshift(fftpkt.fft(self.curr_sig))
        self.curr_freqs = fftpkt.fftshift(fftpkt.fftfreq(len(self.curr_sig), 1.0 / self.sample_rate))

        self.curr_time = (float(peaks[ind]) - filter_size) / self.sample_rate
        self.time_samples_vect.append(self.curr_time)
        if plotting:
            self.PlotPartSig(str(self.curr_time))
        self.prev_time = self.curr_time

class Speaker(object):
    def __init__(self):
        self.id = 555
        self.x = 0
        self.y = 0
        self.z = 0
        self.chirp = []
        self.params = {}
        self.proccessed_signal = SignalHandler()
        self.unfiltered_signal = {}
        self.corr_peaks_1 = []

    def CutSigByPeaks(self, ind, peaks, filter_size):
        if ind >= (len(peaks) - 1):
            self.curr_sig = self.proccessed_signal.filtered_signal[peaks[ind] - filter_size:]
        else:
            final = peaks[ind + 1] - 500
            self.curr_sig = self.proccessed_signal.filtered_signal[peaks[ind] - filter_size: final]

        # calculate FFT
        self.curr_sig_fft = fftpkt.fftshift(fftpkt.fft(self.curr_sig))
        self.curr_freqs = fftpkt.fftshift(fftpkt.fftfreq(len(self.curr_sig), 1.0 / self.sample_rate))

        self.curr_time = (float(peaks[ind]) - filter_size)/self.proccessed_signal.Fs
        if plotting:
            self.PlotPartSig("speaker {0} , current time:{1}".format(self.id, self.curr_time))
        self.prev_time = self.curr_time

    def PlotPartSig(self,title):
        self.part_rec_time = np.linspace(0, float(len(self.curr_sig)) / self.proccessed_signal.Fs, num=len(self.curr_sig))

        plt.figure()
        plt.title('Signal Wave: ' + title )
        plt.plot(self.part_rec_time, self.curr_sig)
        plt.show()

    def BuildChirp(self, freqs_dict, Fs, ch_time, mode, matlab_path):
        # tt = np.linspace(0,0.001, Fs*0.001)
        tt = np.linspace(0, ch_time, int(Fs*ch_time))
        matlab_chirps = io.loadmat(matlab_path)
        chirps = matlab_chirps['allchirp']
        self.unfiltered_signal = {'Fs': Fs, 'low_freq': freqs_dict[str(self.id)][1],
                                  'high_freq': freqs_dict[str(self.id)][0], 'chirp_time': ch_time}

        if mode == 3:
            tmpchirp = signal.chirp(tt, self.unfiltered_signal['low_freq'],
                                    self.unfiltered_signal['chirp_time'],
                                    self.unfiltered_signal['high_freq'])
            self.unfiltered_signal['normalized_vect'] = np.linspace(1, 0.02, num=len(tmpchirp))
            sig1 = np.multiply(self.unfiltered_signal['normalized_vect'], tmpchirp), chirps[:, self.id - 1]
        elif mode == 2:
            ch1 = signal.chirp(tt[0:int(Fs*ch_time/2)], self.unfiltered_signal['low_freq'],
                                    self.unfiltered_signal['chirp_time']/2,
                                    self.unfiltered_signal['high_freq'])
            ch2 = signal.chirp(tt, self.unfiltered_signal['high_freq'],
                                    self.unfiltered_signal['chirp_time']/2,
                                    self.unfiltered_signal['low_freq'])
            tmpchirp = np.concatenate((ch1,ch2))
            self.unfiltered_signal['normalized_vect'] = signal.windows.hamming(len(tmpchirp))
            sig1 = np.multiply(tmpchirp,self.unfiltered_signal['normalized_vect']), chirps[:, self.id - 1]
        elif mode == 1:
            sig1 = signal.chirp(tt, self.unfiltered_signal['low_freq'],
                                    self.unfiltered_signal['chirp_time'],
                                    self.unfiltered_signal['high_freq']), chirps[:, self.id - 1]
        else:
            print "wrong input"

        return sig1

    # if we want to use the waveform builder

    '''
    def WaveformBuilder(self):
        if (('type' not in self.params) | ('Fs' not in self.params) | ('f0' not in self.params) | ('f1' not in self.params) | ('t1' not in self.params) | ('nrm' not in self.params) | ('mode' not in self.params)):
            raise ValueError
        else:
            tt = np.linspace(0, 0.001, self.params['Fs'] * 0.001)
            if self.params['type'] == 'builtin':
                tmp_chirp = signal.chirp(tt, self.params['f0'], self.params['t1'], self.params['f1'],method=self.params['mode'])
            elif self.params['type'] == 'unique':
                k = (self.params['f1'] - self.params['f0']) / (self.params['t1'])
                f = self.params['f0'] + k * tt + cos(tt)
                tmp_chirp = sin(2 * pi * f * tt)
            if self.params['nrm']:
                nrm = np.linspace(1, 0.02, num=len(tmp_chirp))
                final_sig = np.multiply(nrm, tmp_chirp)
            else:
                final_sig = tmp_chirp

            return final_sig

    def Define_ID_WFB(self,my_id):
        self.id = my_id
        self.params['type'] = 'builtin'  # { builtin , unique }
        self.params['t1'] = 0.001
        self.params['Fs'] = 88200
        self.params['mode'] = 'quadratic'  # { logarithmic, linear, hyperbolic , quadratic }
        self.params['nrm'] = False
        if self.id == 1:
            self.params['f0'] = 30000
            self.params['f1'] = 27500
        elif self.id == 2:
            self.params['f0'] = 15000
            self.params['f1'] = 12500
        elif self.id == 3:
            self.params['f0'] = 20000
            self.params['f1'] = 17500
        elif self.id == 4:
            self.params['f0'] = 25000
            self.params['f1'] = 22500

        self.chirp = self.WaveformBuilder()

        matlabchirps = io.loadmat('inputs/all_chirp.mat')
        chirps = matlabchirps['allchirp']
        self.matlab_chirp = chirps[:,self.id]
        
    '''

    def Define_ID(self,my_id, freqs_dict, Fs, ch_time, matlab_path):
        '''
        define all params , generate signals for correlation, load matlab TX signals ,
        filters record following to chirp parameters
        :param my_id: speaker ID
        :param freqs_dict: all expected chirps frequencies
        :param Fs: recording sample_rate.
        :return:
        '''
        self.id = my_id
        self.chirp, self.matlab_chirp = self.BuildChirp(freqs_dict, Fs, ch_time, 1, matlab_path)

    def DefineLocation(self,op,x=0.0,y=0.0,z=0.0):
        if op == 'csv':
            return
        elif op == 'manual':
            self.x = float(raw_input("x:").strip())
            self.y = float(raw_input("y:").strip())
            self.z = float(raw_input("z:").strip())
        else:
            self.x = x
            self.y = y
            self.z = z



if __name__ == '__main__':

    sp_list = [Speaker(), Speaker(), Speaker(), Speaker()]
    # if TOA already exist
    TOA_path = ''    # '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/output/toa_record_1542718040.csv'
    speakers_frequencies = {'1': [20000, 27000],
                            '2': [27000, 34000],
                            '3': [34000, 41000],
                            '4': [42000, 49000]}
    speakers_locations_d = {'1': [0.0, 0.0, 2.30],
                            '2': [3.9, 0.0, 2.25],
                            '3': [0.65, 3.8, 2.23],
                            '4': [3.9, 4.04, 1.52]}
    room_sizes = {'x': 3.9,
                  'y': 4.04,
                  'z': 2.4}

    chirp_time = 0.005
    filter_size = 1001
    matlab_path = '../inputs/chirp_blackman_harris_5ms_Fs500000.mat'
    record_path = '../inputs/blackmanharris5ms/4.WAV'

    # algorithm = int(raw_input("choose algorithm:\n\t(1) Chan Algorithm\n\t"
    #                           "(2) Taylor Algorithm\n\t(3) Room LUT Algorithm\n\t"
    #                           "(4) All Algorithms").strip())
    algorithm = 1
    use_averaging_before_calculation = False
    time_factor = 2
    avg_group_size = 5
    resolution = 0.05           # [m]
    utils_obj = UTILS()
    record = recwav()
    if not (os.path.isfile(TOA_path) & TOA_path.endswith('.csv')):
        record.change_path(os.path.abspath(record_path), 'in')
        record.PlotSignal('blackmanharris5ms')
        record.PlotFFT(record.path)
        # record.Spectogram()

    for i in range(len(sp_list)):
        sp_list[i].Define_ID(i + 1, speakers_frequencies, record.sample_rate, chirp_time, matlab_path)
        sp_list[i].DefineLocation('s',
                                  speakers_locations_d[str(sp_list[i].id)][0],
                                  speakers_locations_d[str(sp_list[i].id)][1],
                                  speakers_locations_d[str(sp_list[i].id)][2])

    if not (os.path.isfile(TOA_path) & TOA_path.endswith('.csv')):
        # define speaker parameters and filtering signals according to speakers frequencies
        for i in range(len(sp_list)):
            # resampling(sp_list[i], record.sample_rate)
            sp_list[i].unfiltered_signal['signal'] = record.signal
            sp_list[i].unfiltered_signal['time_vect'] = record.rec_time
            sp_list[i].proccessed_signal.defineParams(sp_list[i].unfiltered_signal)
            sp_list[i].proccessed_signal.BPF(filter_size,plotting=plotting)
            sp_list[i].peaks, sp_list[i].peaks_height = signal.find_peaks(
                sp_list[i].proccessed_signal.filtered_signal,
                height=int(0.6*(max(sp_list[i].proccessed_signal.filtered_signal))),
                distance=50000)
            sp_list[i].correl1, sp_list[i].correl2 = utils_obj.CorrSpeakerSig(sp_list[i])
            sp_list[i].corr_peaks, sp_list[i].corr_peaks_height = signal.find_peaks(
                sp_list[i].correl1,
                height=int(0.6*(max(sp_list[i].correl1))),
                distance=50000)
            sp_list[i].corr_peaks_test, sp_list[i].corr_peaks_height_test = signal.find_peaks(
                sp_list[i].correl1,
                height=48000,
                distance=500)

        for i in range(len(sp_list[0].corr_peaks)):
            min_peak = min([sp_list[j].corr_peaks[i] for j in range(len(sp_list))])
            # first_sp = np.argmin(ith_peaks)

            for sp in sp_list:
                a = sp.corr_peaks_test - min_peak
                for k in range(len(a)):
                    if a[k] < 0:
                        a[k] = 10000000000
                loac = np.argmin(a)
                # sp.corr_peaks_1.append(sp.corr_peaks_test[loac])
                sp.corr_peaks[i] = sp.corr_peaks_test[loac]


        if plotting:
            print "plotting"
            plt.figure()
            count = 0
            for sp in sp_list:
                count += 1
                plt.subplot(2, 2, count)
                plt.plot(sp.unfiltered_signal['time_vect'], sp.unfiltered_signal['signal'])
                plt.plot(sp.unfiltered_signal['time_vect'],
                         sp.proccessed_signal.filtered_signal[:len(sp.unfiltered_signal['signal'])])
                # plt.plot(sp.peaks, sp.proccessed_signal.filtered_signal[sp.peaks], "x")
                # plt.plot(np.zeros_like(sp.proccessed_signal.filtered_signal), "--", color="gray")
                plt.grid()
            plt.show()
            # count = 0
            # fftd = {}
            # freqsd = {}
            # plt.figure(790)
            # for sp in sp_list:
            #     count += 1
            #     fftd[str(sp.id)] = fftpkt.fftshift(fftpkt.fft(sp.proccessed_signal.filtered_signal))
            #     freqsd[str(sp.id)] = fftpkt.fftshift(fftpkt.fftfreq(len(sp.proccessed_signal.filtered_signal), 1.0 / sp.proccessed_signal.Fs))
            #     plt.subplot(2, 2, count)
            #     plt.plot(freqsd[str(sp.id)], abs(fftd[str(sp.id)]))
            # plt.show()

        with open(record.toa_csv_path, 'wb') as fout:
            writer = csv.DictWriter(fout, fieldnames=["toa_sp_1", "toa_sp_2", "toa_sp_3", "toa_sp_4"])
            writer.writeheader()
            iterr = 1
            peaks, _ = signal.find_peaks(record.signal, height=int(0.7*(max(record.signal))), distance=10000)
            if plotting:
                plt.plot(record.signal)
                plt.plot(peaks, record.signal[peaks], "x")
                plt.plot(np.zeros_like(record.signal), "--", color="gray")
                plt.show()

            # delay = (filter_size - 1) / 2 + int(len(sp_list[0].chirp) / 2)
            _sum = 0
            for sp in sp_list:
                # sp.peaks_after_delay = sp.corr_peaks - delay
                sp.peaks_after_delay = sp.corr_peaks
                sp.peaks_time_stamps = sp.peaks_after_delay / float(sp.proccessed_signal.Fs)
                _sum += sp.peaks_time_stamps
            timestamps = _sum / float(len(sp_list)) - 5 * (10 ** -3)


            # for i in range(len(peaks)):
            #     # record.CutCurrSig()
            #     # record.CutSigByPeaks(i, peaks, filter_size)
            #     # for sp in sp_list:
            #     #     sp.CutSigByPeaks(i, peaks, filter_size)

            # utils_obj.CorrWith4(sp_list)

            corr_time_vect = np.linspace(0, float(len(sp_list[0].correl1)) / sp_list[0].proccessed_signal.Fs,
                                         num=len(sp_list[0].correl1))

            if plotting:
                plt.figure()
                count = 0
                for sp in sp_list:
                    count += 1
                    plt.subplot(2, 2, count)
                    plt.plot(sp.unfiltered_signal['time_vect'], sp.unfiltered_signal['signal'])
                    plt.plot(sp.proccessed_signal.record_time_with_filter_delay,
                             sp.proccessed_signal.filtered_signal)
                    plt.plot(corr_time_vect,sp.correl1)
                    plt.grid()
                    for i in range(len(sp.peaks)):
                        print ("*" * 50 + "\n") * 3
                        # print "real_peak = {0}".format(peaks[i])
                        print "real_peak_filtered = {0}".format(sp.peaks[i])
                        print "corr peak = {0}".format(sp.corr_peaks[i])
                        print "corr_peak - real_peak_filtered = {0}".format(sp.corr_peaks[i] - sp.peaks[i])
                        print ("*" * 50 + "\n") * 3
                plt.show()

                plt.figure()
                for sp in sp_list:
                    plt.plot(corr_time_vect, sp.correl1)
                plt.legend(["sp1", "sp2", "sp3", "sp4"])
                plt.show()
                plt.figure()
                # for sp in sp_list:
                #     plt.plot(sp.correl1)
                #     # plt.plot(sp.correl2)
                # plt.legend(["sp1", "sp2", "sp3", "sp4"])
                # # plt.legend(["sp1", "sp1-no filt", "sp2", "sp2-no filt", "sp3", "sp3-no filt", "sp4", "sp4-no filt"])
                # plt.show()


            for i in range(len(sp_list[0].peaks_time_stamps)):
                writer.writerow({"toa_sp_1": sp_list[0].peaks_time_stamps[i],
                                 "toa_sp_2": sp_list[1].peaks_time_stamps[i],
                                 "toa_sp_3": sp_list[2].peaks_time_stamps[i],
                                 "toa_sp_4": sp_list[3].peaks_time_stamps[i]})

        rec_dict = utils_obj.csv2dict(record.toa_csv_path, {"toa_sp_1": [],
                                                            "toa_sp_2": [],
                                                            "toa_sp_3": [],
                                                            "toa_sp_4": []})

    else:
        record.toa_csv_path = TOA_path
        rec_dict = utils_obj.csv2dict(record.toa_csv_path, {"toa_sp_1": [],
                                                            "toa_sp_2": [],
                                                            "toa_sp_3": [],
                                                            "toa_sp_4": []})

        timestamps = [sum([rec_dict["toa_sp_" + str(i+1)][j]/len(rec_dict) - (5 * 10**(-3))
                           for i in range(len(rec_dict))])
                      for j in range(len(rec_dict["toa_sp_1"]))]





    print colored("Finished parsing wav file","green")

    # throw outliers according to 3D Gaussian model
    if use_averaging_before_calculation:    # throw outliers only by tdoa_from_sp_1 samples calculation
        TDOA_for_outliers = utils_obj.CreateTDOAlistBySPx(1,rec_dict)
        utils_obj.ScatterPlot3D(TDOA_for_outliers['tdoa_sp_2'], TDOA_for_outliers['tdoa_sp_3'], TDOA_for_outliers['tdoa_sp_4'],
                                'TDOA results from sp1', ['TDOA_21 [sec]', 'TDOA_31 [sec]', 'TDOA_41 [sec]'],
                                [(-0.001,-0.003), (-0.003,-0.005), (-0.005,-0.007)])

        rec_dict, avg_list = utils_obj.Throw_Outliers(rec_dict, TDOA_for_outliers, avg_group_size)
    else:
        avg_list = [1] * len(rec_dict["toa_sp_1"])



    print colored("generate sp2mic list", "blue")
    sp2mic = [rec_dict['toa_sp_1'], rec_dict['toa_sp_2'], rec_dict['toa_sp_3'], rec_dict['toa_sp_4']]
    sp_location = utils_obj.buildSpeakersLocationMatrix(sp_list)

    if algorithm == 1:
        chan = algos.ChanAlgo(avg_dim=avg_group_size,
                              use_avg=use_averaging_before_calculation,
                              )
        location_list = chan.chan_main(sp2mic, sp_location, timestamps, avg_list)

    elif algorithm == 2:
        taylor_obj = algos.TaylorLS()
        print "TBD"
        results_dict = {}
    elif algorithm == 3:
        LUT_obj = algos.RoomMatrix(room_sizes['x'],
                                   room_sizes['y'],
                                   room_sizes['z'],
                                   sp_list,
                                   res=0.05,
                                   avg_dim=avg_group_size)
        location_list = LUT_obj.RoomMatMain(sp2mic, timestamps, avg_list,
                                            room_shape='square', use_avg=use_averaging_before_calculation)
    elif algorithm == 4:
        chan = algos.ChanAlgo()
        LUT_obj = algos.RoomMatrix()
        taylor_obj = algos.TaylorLS()

        location_list_chan = chan.chan_main(sp2mic, sp_location)
        location_list_lut = LUT_obj.RoomMatMain()
        location_list_taylor = []  # TBD
        path_chan = os.path.abspath('../output') + '/chan_locations_results_'+str(int(time.time()))+'.csv'
        path_lut = os.path.abspath('../output') + '/room_matrix_locations_results_' + str(int(time.time())) + '.csv'
        # path_taylor = os.path.abspath('../output') + '/taylor_ls_locations_results_' + str(int(time.time())) + '.csv'
        res_dict_chan = utils_obj.res2csv(record.time_samples_vect, location_list_chan,path_chan)
        res_dict_LUT = utils_obj.res2csv(record.time_samples_vect, location_list_lut, path_lut)
        # res_dict_taylor = utils_obj.res2csv(record.time_samples_vect, location_list_taylor, path_taylor)  # TBD


    if not use_averaging_before_calculation:
        from statistics import mean

        if algorithm == 1:
            location_list = chan.FilterRoom(location_list,[(0, room_sizes['x']),
                                                           (0, room_sizes['y']),
                                                           (0, room_sizes['z'])])

        number_of_groups = int(len(location_list)/(record.total_time / time_factor))
        num_of_elements = int(len(location_list)/number_of_groups)

        avg_time = []
        avg_error = []
        avg_location = [[],[],[]]
        for k in range(number_of_groups):
            v = [[], [], []]
            locloc = location_list[k*num_of_elements: k*num_of_elements + num_of_elements]

            for i in range(3):
                v[i] = [loc_[1][i] for loc_ in locloc]

            outliers_list = utils_obj.Find_Outliers(v)

            v_new = [[],[],[]]
            tmp_time = []
            tmp_error = []
            for k in range(len(outliers_list)):
                if outliers_list[k] == False:
                    for i in range(3):
                        v_new[i].append(v[i][k])
                    tmp_time.append(locloc[k][0])
                    tmp_error.append(locloc[k][2])
            avg_time.append(mean(tmp_time))
            if algorithm == 3: # LUT matrix
                avg_error.append(mean(tmp_error))
            else:
                avg_error.append(0)
            for i in range(3):
                avg_location[i].append(mean(v_new[i]))

        location_list = [[avg_time[i], [avg_location[0][i],avg_location[1][i],avg_location[2][i]], avg_error[i]]
                         for i in range(len(avg_time))]

    if algorithm < 4:
        results_dict = utils_obj.res2csv(location_list,record.results_path)

    res2print = df(results_dict)
    print colored(res2print,'green')


    # need to add calculate error from expected value (Euclidean Distance)

    utils_obj.ScatterPlot3D(results_dict['X [m]'],results_dict['Y [m]'], results_dict['Z [m]'],
                            'Algrithm localization decision - location of microphone', ['X[m]', 'Y[m]', 'Z[m]'],
                            [(0, room_sizes['x']), (0, room_sizes['y']), (0, room_sizes['z'])])


