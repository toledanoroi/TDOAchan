# project: TDOA based ultrasonic sound localization system
# authors: Roi Toledano & Yarden Avraham
# lab    : Bats lab
# guide  : PhD Yossi Yovel

import csv
import wave
import numpy as np
from pandas import DataFrame as df
import time
import os
from scipy import signal, io
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
        plt.show(block=False)

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
        plt.show(block=False)


class recwav(object):
    def __init__(self, point_name):
        self.path = '../inputs/blackmanharris5ms/1.wav'
        self.sample_rate = 96000
        self.toa_csv_path = '../output/' + point_name + '.csv'
        self.results_path = '../output/locations_results_' + str(int(time.time()))+'.csv'
        self.record_sig = []
        self.num_of_channels = 1
        self.prev_time = 0
        self.curr_time = 0
        self.parttime = 0.49
        self.total_time = 0
        self.time_samples_vect = []

    def PlotSignal(self, title):
        '''
        reads signal wav file to the python as ndarray, update record parameters according to that signal.
        and plot it.
        :param title: the title of the plot.
        :return:
        '''
        self.sig_file_pointer = wave.open(self.path, 'r')
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
            plt.figure()
            plt.title('Signal Wave: ' + title )
            plt.plot(self.rec_time, self.signal)
            plt.show(block=False)

    def PlotFFT(self,title):
        '''
        calculate the FFT and psd of the record and plot it.
        :param title: the title of the plot.
        '''
        self.fft = fftpkt.fft(self.signal)
        self.psd = np.abs(self.fft) ** 2
        self.freqs = fftpkt.fftshift(fftpkt.fftfreq(len(self.signal),1.0/self.sample_rate))
        self.fft = fftpkt.fftshift(self.fft)

        if plotting:
            plt.figure()
            plt.title('Signal Wave FFT: ' + title)
            plt.plot(self.freqs, abs(self.fft))
            plt.grid()
            plt.show(block=False)

    def change_path(self, path, mode='in'):
        '''
        change one of the record paths
        :param path: string of existing file
        :param mode: in for record WAV file else for toa csv path name
        '''

        if mode == 'in':
            if os.path.isfile(path):
                self.path = path
            else:
                print "the path isn't exist"
                raise ValueError()
        else:
            self.toa_csv_path = path

    def Spectogram(self):
        '''
        Calculate the spectogram of the signal and plot it.
        :return:
        '''
        f, t, Sxx = signal.spectrogram(self.signal,self.sample_rate, nperseg=64)
        if plotting:
            plt.pcolormesh(t, f, (10 * np.log10(Sxx)), cmap='inferno')
            plt.colorbar()
            plt.ylabel('Frequency [Hz]')
            plt.xlabel('Time [sec]')
            plt.show(block=False)
        return f, t, Sxx



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

    def BuildChirp(self, freqs_dict, Fs, ch_time, mode, matlab_path):
        '''
        generate chirp in python according to signal parameters
        and load matlab expected signal to python as ndarray
        :param freqs_dict: frequencies of the chirp
        :param Fs: sample rate
        :param ch_time: chirp duration time
        :param mode: which type of chirp
        :param matlab_path: matlab expected signal .mat file path
        :return: python generated expected signal and matlab generated expected signal as tuple
        '''
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

    def Define_ID(self, my_id, freqs_dict, Fs, ch_time, matlab_path, signal_mode):
        '''
        define all params , generate signals for correlation, load matlab TX signals ,
        filters record following to chirp parameters
        :param my_id: speaker ID
        :param freqs_dict: all expected chirps frequencies
        :param Fs: recording sample_rate.
        :param matlab_path: path of the expected signal, generated by the transmitting unit (matlab)
        :param signal_mode: type of signal
        :return:
        '''
        self.id = my_id
        self.chirp, self.matlab_chirp = self.BuildChirp(freqs_dict, Fs, ch_time, signal_mode, matlab_path)

    def DefineLocation(self, op, x=0.0, y=0.0, z=0.0):
        '''
        define the location of the speaker
        :param op: how to define it ? {'csv' , 'manual' or from code }
        :param x: location in x axis
        default: x = 0.0
        :param y: location in y axis
        default: y = 0.0
        :param z: location in z axis
        default: z = 0.0
        :return:
        '''
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


def RxMain(params):
    '''
    The main funtion of the Receiver (microphone) signal processing
    :param params: parameters of the test
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
            :param TOA_path - if toa samples already exists: path of the relevant TOA samples. else = ''
            :type string , should br a real path of csv file
            ---------------------------------------------------------------------------------------------
            :param matlab_path - path of the expected signals that generated by the transmitting unit
            (matlab)
            :type string , should br a real path of mat file with variable 'allchirp'
            ---------------------------------------------------------------------------------------------
            :param record_path - path of the recorded signal to examine.
            :type string , should br a real path of .WAV file
            ---------------------------------------------------------------------------------------------
            :param speakers_frequencies the recorded chirp frequencies that defined in the Transmitting
            unit.
            stracture: {'1': [start_freq, stop_freq],'2': [start_freq, stop_freq],...}
            :type  dictionary
            ---------------------------------------------------------------------------------------------
            :param speakers_locations_d - locations dictionary for all the speakers
            stracture: {'1': [x, y, z],'2': [x, y, z],...}
            :type dictionary of int lists
            ---------------------------------------------------------------------------------------------
            :param room_sizes - dictionary of the room limits with every axis
            stracture: {'x':<maxX>,'y':<maxY>,'z':<maxZ>}
            :type
            ---------------------------------------------------------------------------------------------
            :param one_shot - define if to take one record and parsed it or take toacsv.csv with ready toa paths
            :type boolean , true for .wav file False for .csv files
            ---------------------------------------------------------------------------------------------
            :param res_iteration - true if want to examine the matrix resolution feature
            else define a specific resolution
            :type boolean
            ---------------------------------------------------------------------------------------------
            :param only_toa - define if to run the Rxmain with algorithm calculation or only TOA measurements
            True -> only TOA measurements , False -> with algorithm calculation.
            :type boolean
            ---------------------------------------------------------------------------------------------
            :param chirp_time - the duration of the chirp signal generated by matlab
            :type double
            ---------------------------------------------------------------------------------------------
            :param filter_size - number of taps in the speakers Band Pass Filters (BPF)
            :type int , must be odd number for real BPF (FIR design)
            ---------------------------------------------------------------------------------------------
            :param signal_mode
            :type
            ---------------------------------------------------------------------------------------------
            :param expected_points - an array of expected (x,y,z) defined by human measurements.
            :type
            ---------------------------------------------------------------------------------------------
            :param algorithm - which algorithm to run in the test, choose from the algorithms dictionary:
            algorithm_d = {'chan': 1,'taylor': 2,'room': 3,'both': 4}
            :type string
            ---------------------------------------------------------------------------------------------
            :param resolution - the LUT room matrix algorithm resolution
            :type float
            ---------------------------------------------------------------------------------------------
            :param use_averaging_before_calculation - define where to search for outliers and average results,
            before calculate distance or after.
            :type boolean , True -> before , False -> after
            ---------------------------------------------------------------------------------------------
            :param time_factor - if decide after , define the time period to averaging and throw outliers
            :type float
            ---------------------------------------------------------------------------------------------
            :param avg_group_size - if decide before , define how much samples to averaging
            and throw outliers each iteration.
            :type int
            ---------------------------------------------------------------------------------------------
            :param room3D - array 3D of all edges points in the room
            :type array 3D
            ---------------------------------------------------------------------------------------------
            :param triangle3D - if there is edges that doesn't participate as edge in convex hull,
            define array 3D of the these points and theres neighbors
            :type array 3D
            ---------------------------------------------------------------------------------------------
    :return: average errors of the algorithm
    '''
    # validity check
    if ((params['number_of_speakers'] != len(params['speakers_frequencies'])) |
            (params['number_of_speakers'] != len(params['speakers_locations_d']))):
        print colored("Wrong Input Error \n\tnumber_of_speakers = {0}"
                      "\n\tlen(speakers_frequencies) = {1}"
                      "\n\tlen(speakers_locations_d) = {2}".format(params['number_of_speakers'],
                                                                   len(params['speakers_frequencies']),
                                                                   len(params['speakers_locations_d'])), "red")
    sp_list = []
    for j in xrange(params['number_of_speakers']):
        sp_list.append(Speaker())

    utils_obj = UTILS()
    record = recwav(params['point_name'])
    if params['unique_path'] != 'no_path':
        record.results_path = params['unique_path']

    if params['ToAs_file'] != 'no_path':
        record.toa_csv_path = params['ToAs_file']


    # define room convex
    hull, non_hull, hull2d, non_hull2d = utils_obj.DefineRoom(params['room3D'],
                                                             params['triangle3D'],
                                                             shapewithnonedgepoint=True,
                                                             plotting=False
                                                             )
    if not (os.path.isfile(params['TOA_path']) & params['TOA_path'].endswith('.csv')):
        record.change_path(os.path.abspath(params['record_path']), 'in')
        record.PlotSignal('blackmanharris5ms')
        record.PlotFFT(record.path)
        spectf, spectt, spectsxx = record.Spectogram()
    for i in range(len(sp_list)):
        sp_list[i].Define_ID(i + 1,
                             params['speakers_frequencies'],
                             record.sample_rate,
                             params['chirp_time'],
                             params['matlab_path'],
                             params['signal_mode']
                             )
        sp_list[i].DefineLocation('s',
                                  params['speakers_locations_d'][str(sp_list[i].id)][0],
                                  params['speakers_locations_d'][str(sp_list[i].id)][1],
                                  params['speakers_locations_d'][str(sp_list[i].id)][2]
                                  )
    if not (os.path.isfile(params['TOA_path']) & params['TOA_path'].endswith('.csv')):
        # define speaker parameters and filtering signals according to speakers frequencies
        peaks, _ = signal.find_peaks(record.signal, distance=int(0.8 * record.sample_rate / 2))
        if plotting:
            plt.plot(record.signal)
            plt.plot(peaks, record.signal[peaks], "x")
            plt.plot(np.zeros_like(record.signal), "--", color="gray")
            plt.show(block=False)
        #     here is where we handle the peaks finder
        # --------------------------------------------------------
        # --------------------------------------------------------
        # --------------------------------------------------------
            # --------------------------------------------------------
            # --------------------------------------------------------
            # --------------------------------------------------------
        for i in range(len(sp_list)):
            # resampling(sp_list[i], record.sample_rate)
            sp_list[i].unfiltered_signal['signal'] = record.signal
            sp_list[i].unfiltered_signal['time_vect'] = record.rec_time
            sp_list[i].proccessed_signal.defineParams(sp_list[i].unfiltered_signal)
            sp_list[i].proccessed_signal.BPF(params['filter_size'], plotting=plotting)
            sp_list[i].peaks, sp_list[i].peaks_height = signal.find_peaks(
                sp_list[i].proccessed_signal.filtered_signal,
                # height=int(0.4 * (max(sp_list[i].proccessed_signal.filtered_signal))),
                distance=int(0.9 * sp_list[i].proccessed_signal.Fs / 2)
            )
            sp_list[i].correl1, sp_list[i].correl2 = utils_obj.CorrSpeakerSig(sp_list[i])
            sp_list[i].corr_peaks, sp_list[i].corr_peaks_height = signal.find_peaks(
                sp_list[i].correl1,
                # height=int(0.4 * (max(sp_list[i].correl1))),
                distance=int(0.9 * sp_list[i].proccessed_signal.Fs / 2)    #approximatly 0.45 seconds
            )

            sp_list[i].tot_corr, sp_list[i].tot_corr_height = signal.find_peaks(
                sp_list[i].correl1,
                height=4000) # noise correlation height threshold
            sp_list[i].half = np.median(sp_list[i].tot_corr_height['peak_heights'])

            sp_list[i].corr_peaks_test, sp_list[i].corr_peaks_height_test = signal.find_peaks(
                sp_list[i].correl1,
                height=sp_list[i].half,
                distance=400) # distance between consecutive peaks

            # --------------------------------------------------------
            # --------------------------------------------------------
            # --------------------------------------------------------
            # --------------------------------------------------------
            # --------------------------------------------------------
            # --------------------------------------------------------
        # align corr peaks lists to the same size
        tmp = min([len(sp.corr_peaks) for sp in sp_list])
        for sp in sp_list:
            if len(sp.corr_peaks) > tmp:
                sp.corr_peaks = sp.corr_peaks[0:tmp]

        if plotting:
            plt.figure()
            plt.plot(sp_list[0].correl1)
            plt.plot(sp_list[0].corr_peaks, sp_list[0].correl1[sp_list[0].corr_peaks], "x")
            plt.plot(np.zeros_like(sp_list[0].correl1), "--", color="gray")
            plt.show(block=False)
            plt.figure()
            plt.plot(sp_list[1].correl1)
            plt.plot(sp_list[1].corr_peaks, sp_list[1].correl1[sp_list[1].corr_peaks], "x")
            plt.plot(np.zeros_like(sp_list[1].correl1), "--", color="gray")
            plt.show(block=False)
            plt.figure()
            plt.plot(sp_list[2].correl1)
            plt.plot(sp_list[2].corr_peaks, sp_list[2].correl1[sp_list[2].corr_peaks], "x")
            plt.plot(np.zeros_like(sp_list[1].correl1), "--", color="gray")
            plt.show(block=False)
            plt.figure()
            plt.plot(sp_list[3].correl1)
            plt.plot(sp_list[3].corr_peaks, sp_list[3].correl1[sp_list[3].corr_peaks], "x")
            plt.plot(np.zeros_like(sp_list[1].correl1), "--", color="gray")
            plt.show(block=False)

        for i in range(len(sp_list[0].corr_peaks)):
            min_peak = min([sp_list[j].corr_peaks[i] for j in range(len(sp_list))])

            for sp in sp_list:
                a = sp.corr_peaks_test - min_peak
                for k in range(len(a)):
                    if a[k] < -6000:
                        a[k] = 10000000000
                loac = np.argmin(a)
                sp.corr_peaks[i] = sp.corr_peaks_test[loac]

        if plotting:
            print "plotting"
            plt.figure()
            count = 0
            for sp in sp_list:
                count += 1
                plt.subplot(2, 2, count)
                plt.title("Speaker " + str(sp.id) + " - signal vs filtered signal")
                plt.plot(sp.unfiltered_signal['time_vect'], sp.unfiltered_signal['signal'])
                plt.plot(sp.unfiltered_signal['time_vect'],
                         sp.proccessed_signal.filtered_signal[:len(sp.unfiltered_signal['signal'])])
                plt.grid()
            plt.show(block=False)

        toa_texts = ["toa_sp_" + str(sp.id) for sp in sp_list]
        with open(record.toa_csv_path, 'wb') as fout:
            writer = csv.DictWriter(fout, fieldnames=toa_texts)
            writer.writeheader()
            iterr = 1
            peaks, _ = signal.find_peaks(record.signal, height=int(0.7 * (max(record.signal))), distance=10000)
            if plotting:
                plt.plot(record.signal)
                plt.plot(peaks, record.signal[peaks], "x")
                plt.plot(np.zeros_like(record.signal), "--", color="gray")
                plt.show(block=False)

            # delay = (params['filter_size'] - 1) / 2 + int(len(sp_list[0].chirp) / 2)
            _sum = 0
            for sp in sp_list:
                # sp.peaks_after_delay = sp.corr_peaks - delay
                sp.peaks_after_delay = sp.corr_peaks
                sp.peaks_time_stamps = sp.peaks_after_delay / float(sp.proccessed_signal.Fs)
                _sum += sp.peaks_time_stamps
            timestamps = _sum / float(len(sp_list)) - 5 * (10 ** -3)




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
                    plt.plot(corr_time_vect, sp.correl1)
                    plt.grid()
                    # for i in range(len(sp.peaks)):
                    #     print ("*" * 50 + "\n") * 3
                    #     # print "real_peak = {0}".format(peaks[i])
                    #     print "real_peak_filtered = {0}".format(sp.peaks[i])
                    #     print "corr peak = {0}".format(sp.corr_peaks[i])
                    #     print "corr_peak - real_peak_filtered = {0}".format(sp.corr_peaks[i] - sp.peaks[i])
                    #     print ("*" * 50 + "\n") * 3
                plt.show(block=False)

                plt.figure()
                for sp in sp_list:
                    plt.plot(corr_time_vect, sp.correl1)
                # plt.legend(["sp1", "sp2", "sp3", "sp4"])
                plt.show(block=False)
                #plt.figure()
                # for sp in sp_list:
                #     plt.plot(sp.correl1)
                #     # plt.plot(sp.correl2)
                # plt.legend(["sp1", "sp2", "sp3", "sp4"])
                # # plt.legend(["sp1", "sp1-no filt", "sp2", "sp2-no filt", "sp3", "sp3-no filt", "sp4", "sp4-no filt"])
                # plt.show(block=False)

            for i in range(len(sp_list[0].peaks_time_stamps)):
                dicto = {"toa_sp_" + str(sp.id): sp.peaks_time_stamps[i] for sp in sp_list}
                writer.writerow(dicto)
        dicto1 = {"toa_sp_" + str(sp.id): [] for sp in sp_list}
        rec_dict = utils_obj.csv2dict(record.toa_csv_path, my_dict=dicto1)
    else:
        dicto1 = {"toa_sp_" + str(sp.id): [] for sp in sp_list}
        record.toa_csv_path = params['TOA_path']
        rec_dict = utils_obj.csv2dict(record.toa_csv_path, my_dict=dicto1)
        record.total_time = int(np.average([rec_dict['toa_sp_' + str(sp.id)][-1] for sp in sp_list]))
        timestamps = [sum([rec_dict["toa_sp_" + str(i + 1)][j] / len(rec_dict) - (5 * 10 ** (-3))
                           for i in range(len(rec_dict))])
                      for j in range(len(rec_dict["toa_sp_1"]))]

    print colored("Finished parsing wav file", "green")


    if params['only_toa']:
        print (("*" * 50) + '\n') * 3
        print "\nonly TOA decided , without the " \
              "algorithm calculation\nTOA path : {}\n".format(record.toa_csv_path)
        print (("*" * 50) + '\n') * 3
        return 0, 0

    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------

    # throw outliers according to 3D Gaussian model
    if params['use_averaging_before_calculation']:  # throw outliers only by tdoa_from_sp_1 samples calculation
        print colored("throw outliers according to 3D Gaussian model", "blue")
        TDOA_for_outliers = utils_obj.CreateTDOAlistBySPx(1, rec_dict)
        if plotting:
            if params['number_of_speakers'] == 4 :
                utils_obj.ScatterPlot3D(TDOA_for_outliers['tdoa_sp_2'], TDOA_for_outliers['tdoa_sp_3'],
                                        TDOA_for_outliers['tdoa_sp_4'],
                                        'TDOA results from sp1', ['TDOA_21 [sec]', 'TDOA_31 [sec]', 'TDOA_41 [sec]'],
                                        [(0.01, -0.01), (0.01, -0.01), (0.01, -0.01)])

        rec_dict, avg_list = utils_obj.Throw_Outliers(rec_dict, TDOA_for_outliers, params['avg_group_size'])
    else:
        avg_list = [1] * len(rec_dict["toa_sp_1"])

    print colored("generate sp2mic list", "blue")
    sp2mic = [rec_dict['toa_sp_' + str(sp.id)] for sp in sp_list]

    sp_location = utils_obj.buildSpeakersLocationMatrix(sp_list)

    if params['algorithm'] == 1:
        chan = algos.ChanAlgo(avg_dim=params['avg_group_size'],
                              use_avg=params['use_averaging_before_calculation'],
                              Temprature_meas=params['Temperature']
                              )
        location_list = chan.chan_main(sp2mic, sp_location, timestamps, avg_list)
    elif params['algorithm'] == 2:
        taylor_obj = algos.TaylorLS()
        print "TBD"
        results_dict = {}
    elif params['algorithm'] == 3:
        LUT_obj = algos.RoomMatrix(params['room_sizes']['x'],
                                   params['room_sizes']['y'],
                                   params['room_sizes']['z'],
                                   sp_list,
                                   res=params['resolution'],
                                   avg_dim=params['avg_group_size'],
                                   constant_z=params['constant_z'],
                                   Temprature_meas=params['Temperature']
                                   )
        location_list = LUT_obj.RoomMatMain(sp2mic,
                                            timestamps,
                                            avg_list,
                                            room_shape='square',
                                            use_avg=params['use_averaging_before_calculation'])
    elif params['algorithm'] == 4:
        chan = algos.ChanAlgo()
        LUT_obj = algos.RoomMatrix()
        taylor_obj = algos.TaylorLS()

        location_list_chan = chan.chan_main(sp2mic, sp_location)
        location_list_lut = LUT_obj.RoomMatMain()
        location_list_taylor = []  # TBD
        path_chan = os.path.abspath('../output') + '/chan_locations_results_' + str(int(time.time())) + '.csv'
        path_lut = os.path.abspath('../output') + '/room_matrix_locations_results_' + str(int(time.time())) + '.csv'
        # path_taylor = os.path.abspath('../output') + '/taylor_ls_locations_results_' + str(int(time.time())) + '.csv'
        res_dict_chan = utils_obj.res2csv(record.time_samples_vect, location_list_chan, path_chan)
        res_dict_LUT = utils_obj.res2csv(record.time_samples_vect, location_list_lut, path_lut)
        # res_dict_taylor = utils_obj.res2csv(record.time_samples_vect, location_list_taylor, path_taylor)  # TBD

    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # averaging results according to 3D gaussian models
    if not params['use_averaging_before_calculation']:
        from statistics import mean
        if params['algorithm'] == 1:
            location_list = chan.FilterRoom(location_list, [(0, params['room_sizes']['x']),
                                                            (0, params['room_sizes']['y']),
                                                            (0, params['room_sizes']['z'])])
        if len(location_list) > (float(record.total_time) / params['time_factor']):
            number_of_groups = int((float(record.total_time) / params['time_factor']))
            num_of_elements = int(len(location_list) / number_of_groups)

            avg_time = []
            avg_error = []
            avg_location = [[], [], []]
            for k in range(number_of_groups):
                v = [[], [], []]
                locloc = location_list[k * num_of_elements: k * num_of_elements + num_of_elements]

                for i in range(3):
                    v[i] = [loc_[1][i] for loc_ in locloc]

                outliers_list = utils_obj.Find_Outliers(v)

                v_new = [[], [], []]
                tmp_time = []
                tmp_error = []
                for k in range(len(outliers_list)):
                    if outliers_list[k] == False:
                        for i in range(3):
                            v_new[i].append(v[i][k])
                        tmp_time.append(locloc[k][0])
                        tmp_error.append(locloc[k][2])
                avg_time.append(mean(tmp_time))
                if params['algorithm'] == 3:  # LUT matrix
                    avg_error.append(mean(tmp_error))
                else:
                    avg_error.append(0)
                for i in range(3):
                    avg_location[i].append(mean(v_new[i]))

            location_list = [[avg_time[i], [avg_location[0][i], avg_location[1][i], avg_location[2][i]], avg_error[i]]
                             for i in range(len(avg_time))]

    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    if (os.path.isfile(params['TOA_path']) & params['TOA_path'].endswith('.csv')):
        point_number = int(params['TOA_path'][params['TOA_path'].rfind('/') + 2: params['TOA_path'].rfind('.')])
    elif (os.path.isfile(params['ToAs_file']) & params['ToAs_file'].endswith('.csv')):
        point_number = int(params['ToAs_file'][params['ToAs_file'].rfind('/') + 2: params['ToAs_file'].rfind('.')])
    else:
        point_number = 1

    if params['algorithm'] < 4:
        results_dict = utils_obj.res2csv(location_list,
                                         record.results_path,
                                         1,
                                         params['expected_points'][0][0],
                                         params['expected_points'][0][1],
                                         params['expected_points'][0][2],
                                         point_number)

    res2print = df(results_dict)
    print colored(res2print, 'green')
    # error calculation
    err2d = []
    err3d = []
    for i in range(len(results_dict['X [m]'])):
        err2d.append(utils_obj.CalcErrorFromExpected([results_dict['X [m]'][i],
                                              results_dict['Y [m]'][i]],
                                             params['expected_points2d']))
        err3d.append(utils_obj.CalcErrorFromExpected([results_dict['X [m]'][i],
                                              results_dict['Y [m]'][i],
                                              results_dict['Z [m]'][i]],
                                             params['expected_points']))

    print "For 2D:\nmax error = {0}\nmin error = {1}\naverage error = {2}".format(max(err2d), min(err2d), np.average(err2d))
    print "For 3D:\nmax error = {0}\nmin error = {1}\naverage error = {2}".format(max(err3d), min(err3d), np.average(err3d))


    if params['algorithm'] < 4:
        results_dict = utils_obj.res2csv(location_list,
                                         record.results_path,
                                         np.average(err3d),
                                         params['expected_points'][0][0],
                                         params['expected_points'][0][1],
                                         params['expected_points'][0][2],
                                         point_number)

    res2print = df(results_dict)
    print colored(res2print, 'green')


    if plotting:
        utils_obj.ScatterPlot3D(results_dict['X [m]'],
                                results_dict['Y [m]'],
                                results_dict['Z [m]'],
                                'Algrithm localization decision - location of microphone 3D',
                                ['X[m]', 'Y[m]', 'Z[m]'],
                                [(0, params['room_sizes']['x']),
                                 (0, params['room_sizes']['y']),
                                 (0, params['room_sizes']['z'])],
                                cvx1=hull,
                                cvx2=non_hull,
                                expected=params['expected_points']
                                )

        utils_obj.ScatterPlot2D(results_dict['X [m]'],
                                results_dict['Y [m]'],
                                'Algrithm localization decision - location of microphone 2D',
                                ["X[m]", "Y[m]"],
                                [(-0.2, params['room_sizes']['x']+0.2),
                                 (-0.2, params['room_sizes']['y']+0.2)],
                                cvx1=hull2d,
                                cvx2=non_hull2d,
                                expected=params['expected_points2d']
                                )

    return np.average(err2d), np.average(err3d)


# old stuff

# speaker : if we want to use the waveform builder

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
    # def CutSigByPeaks(self, ind, peaks, filter_size):
    #     if ind >= (len(peaks) - 1):
    #         self.curr_sig = self.proccessed_signal.filtered_signal[peaks[ind] - filter_size:]
    #     else:
    #         final = peaks[ind + 1] - 500
    #         self.curr_sig = self.proccessed_signal.filtered_signal[peaks[ind] - filter_size: final]
    #
    #     # calculate FFT
    #     self.curr_sig_fft = fftpkt.fftshift(fftpkt.fft(self.curr_sig))
    #     self.curr_freqs = fftpkt.fftshift(fftpkt.fftfreq(len(self.curr_sig), 1.0 / self.sample_rate))
    #
    #     self.curr_time = (float(peaks[ind]) - filter_size)/self.proccessed_signal.Fs
    #     if plotting:
    #         self.PlotPartSig("speaker {0} , current time:{1}".format(self.id, self.curr_time))
    #     self.prev_time = self.curr_time
    #
    # def PlotPartSig(self,title):
    #     self.part_rec_time = np.linspace(0, float(len(self.curr_sig)) / self.proccessed_signal.Fs, num=len(self.curr_sig))
    #
    #     plt.figure()
    #     plt.title('Signal Wave: ' + title )
    #     plt.plot(self.part_rec_time, self.curr_sig)
    #     plt.show(block=False)


    # record functions
    # def CutCurrSig(self):
    #     self.total_time = self.rec_time[-1]
    #     self.curr_time += self.parttime
    #     self.time_samples_vect.append(self.curr_time)
    #     if self.curr_time > self.total_time:
    #         print "curr_time exceed total time , update part to total time"
    #         self.curr_time = self.total_time
    #     self.curr_sig = self.signal[int(self.prev_time*self.sample_rate): int(self.curr_time*self.sample_rate)]
    #     self.PlotPartSig(str(self.curr_time))
    #     self.prev_time = self.curr_time
    #     return self.curr_sig
    #
    # def CutSigByPeaks(self, ind, peaks, filter_size):
    #     if ind >= (len(peaks) - 1):
    #         self.curr_sig = self.signal[peaks[ind] - filter_size:]
    #     else:
    #         final = peaks[ind + 1] - 500
    #         self.curr_sig = self.signal[peaks[ind] - filter_size : final]
    #     # calculate FFT
    #     self.curr_sig_fft = fftpkt.fftshift(fftpkt.fft(self.curr_sig))
    #     self.curr_freqs = fftpkt.fftshift(fftpkt.fftfreq(len(self.curr_sig), 1.0 / self.sample_rate))
    #
    #     self.curr_time = (float(peaks[ind]) - filter_size) / self.sample_rate
    #     self.time_samples_vect.append(self.curr_time)
    #     if plotting:
    #         self.PlotPartSig(str(self.curr_time))
    #     self.prev_time = self.curr_time
    # def PlotPartSig(self,title):
    #
    #     # If Stereo
    #     if self.num_of_channels == 2:
    #         print 'Our Application support only mono files'
    #         sys.exit(0)
    #
    #     self.part_rec_time = np.linspace(0, float(len(self.curr_sig)) / self.sample_rate, num=len(self.curr_sig))
    #
    #     plt.figure(12)
    #     plt.title('Signal Wave: ' + title )
    #     plt.plot(self.part_rec_time,self.curr_sig)
    #     plt.show(block=False)
    # def plotOnlyOneCycle(self):
    #     tmpsig = self.signal[int(0.5*self.sample_rate):int(1.3*self.sample_rate)]
    #     tmpfft = fftpkt.fftshift(fftpkt.fft(tmpsig))
    #     tmppsd = np.abs(tmpfft) ** 2
    #     tmpfreqs = fftpkt.fftshift(fftpkt.fftfreq(len(tmpsig), 1.0 / self.sample_rate))
    #     tmprec_time = np.linspace(0, len(tmpsig) / self.sample_rate, num=len(tmpsig))
    #
    #     if plotting:
    #         plt.figure(3)
    #         plt.title('one cycle')
    #         plt.plot(tmprec_time,tmpsig)
    #         plt.show(block=False)
    #
    #         plt.figure(4)
    #         plt.title('one cycle FFT:')
    #         plt.plot(tmpfreqs,abs(tmpfft))
    #         plt.grid()
    #         plt.show(block=False)