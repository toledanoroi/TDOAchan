import csv
import wave
import numpy as np
import time
import os
from scipy import signal
from scipy import io
import scipy.fftpack as fftpkt
from matplotlib import pyplot as plt
import sys
from termcolor import colored
import src.algorithms as algos
from math import cos, sin, pi
from utils import UTILS

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
    freqs  = fftpkt.fftshift(fftpkt.fftfreq(len(speaker.matlab_chirp), 1.0 / 88200))
    fft22  = fftpkt.fftshift(fftpkt.fft(f))
    freqs2 = fftpkt.fftshift(fftpkt.fftfreq(len(f), 1.0 / 88200))
    fft33  = fftpkt.fftshift(fftpkt.fft(speaker.chirp))
    freqs3 = fftpkt.fftshift(fftpkt.fftfreq(len(speaker.chirp), 1.0 / 88200))

    if plotting:
        plt.figure(10)
        plt.plot(freqs, abs(fft1), 'r')
        plt.plot(freqs2, abs(fft22), 'b')
        plt.plot(freqs3, abs(fft33), 'g')
        plt.legend(['matlab 250000', 'resampled', 'python generated'], loc='best')
        plt.show()


class recwav(object):
    def __init__(self):
        self.path = '../inputs/first_test.wav'
        self.sample_rate = 88200
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
            print 'Our Application support only mono files'
            sys.exit(0)

        self.rec_time = np.linspace(0, len(self.signal) / self.sample_rate , num=len(self.signal))
        self.total_time = self.rec_time[-1]

        print "Signal Parameters:\n\tFile Name:\t{0}\n\tSample rate:\t{1}" \
              "\n\tNumber of channels:\t{2}\n\tDuration Time:\t{3}".format(self.path,
                                                                           self.sample_rate,
                                                                           self.num_of_channels,
                                                                           self.rec_time[-1])

        if plotting:
            plt.figure(1)
            plt.title('Signal Wave: ' + title )
            plt.plot(self.rec_time,self.signal)
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

        plt.figure(2)
        plt.title('Signal Wave FFT: ' + title)
        plt.plot(self.freqs,abs(self.fft))
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
        # self.PlotPartSig(str(self.curr_time))
        self.prev_time = self.curr_time
        return self.curr_sig


class Speaker(object):
    def __init__(self):
        self.id = 555
        self.x = 0
        self.y = 0
        self.z = 0
        self.chirp = []
        self.params = {}

    def BuildChirp(self):
        Fs = 88200
        tt = np.linspace(0,0.001,Fs*0.001)
        matlabchirps = io.loadmat('../inputs/all_chirp.mat')
        chirps = matlabchirps['allchirp']
        if self.id == 1:
            tmpchirp = signal.chirp(tt,30000,0.001,27500)
            nrm = np.linspace(1,0.02,num=len(tmpchirp))
            return np.multiply(nrm,tmpchirp) , chirps[:,0]
        elif self.id == 2:
            tmpchirp = signal.chirp(tt,15000,0.001,12500)
            nrm = np.linspace(1,0.02,num=len(tmpchirp))
            return np.multiply(nrm,tmpchirp) , chirps[:,1]
        elif self.id == 3:
            tmpchirp = signal.chirp(tt,20000,0.001,17500)
            nrm = np.linspace(1,0.02,num=len(tmpchirp))
            return np.multiply(nrm,tmpchirp) , chirps[:,2]
        elif self.id == 4:
            tmpchirp = signal.chirp(tt,25000,0.001,22500)
            nrm = np.linspace(1,0.02,num=len(tmpchirp))
            return np.multiply(nrm,tmpchirp) , chirps[:,3]

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

    def Define_ID(self,my_id):
        # define all params , generate signals for correlation, load matlab TX signals
        self.id = my_id
        # self.params['type'] = 'builtin'  # { builtin , unique }
        # self.params['t1'] = 0.001
        # self.params['Fs'] = 88200
        # self.params['mode'] = 'quadratic'  # { logarithmic, linear, hyperbolic , quadratic }
        # self.params['nrm'] = False
        # if self.id == 1:
        #     self.params['f0'] = 30000
        #     self.params['f1'] = 27500
        # elif self.id == 2:
        #     self.params['f0'] = 15000
        #     self.params['f1'] = 12500
        # elif self.id == 3:
        #     self.params['f0'] = 20000
        #     self.params['f1'] = 17500
        # elif self.id == 4:
        #     self.params['f0'] = 25000
        #     self.params['f1'] = 22500
        #
        # self.chirp = self.WaveformBuilder()
        #
        # matlabchirps = io.loadmat('inputs/all_chirp.mat')
        # chirps = matlabchirps['allchirp']
        # self.matlab_chirp = chirps[:,self.id]

        self.chirp, self.matlab_chirp = self.BuildChirp()

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

    def GenerateTX(self):
        # add matlab communicator
        pass



if __name__ == '__main__':

    correl = []
    sp_list = [Speaker(), Speaker(), Speaker(), Speaker()]

    # for item in sp_list:
    #     item.DefineLocation('manual')

    sp_list[0].DefineLocation('s', 0.0, 0.0, 1.0)
    sp_list[1].DefineLocation('s', 0.56346, 3.1346, 1.35)
    sp_list[2].DefineLocation('s', 2.1560, 0.3170, 1.370)
    sp_list[3].DefineLocation('s', 2.10, 3.15730, 1.45)
    algorithm = int(raw_input("choose algorithm:\n\t(1) Chan Algorithm\n\t"
                              "(2) Taylor Algorithm\n\t(3) Room LUT Algorithm\n\t"
                              "(4) All Algorithms").strip())
    utils_obj = UTILS()
    record = recwav()
    record.change_path('../inputs/second_test.wav','in')
    record.PlotSignal('second.wav')
    # record.PlotPartSig('second.wav')
    # record.PlotFFT('second.wav')
    # record.plotOnlyOneCycle()
    # record.Spectogram()

    for i in range(len(sp_list)):
        sp_list[i].Define_ID(i+1)
        resampling(sp_list[i], record.sample_rate)

    with open(record.toa_csv_path, 'wb') as fout:
        writer = csv.DictWriter(fout, fieldnames=["iteration","toa_sp_1", "toa_sp_2", "toa_sp_3",
                                                  "toa_sp_4", "corr_1","corr_2","corr_3","corr_4"])
        writer.writeheader()
        iterr = 1
        while(record.curr_time<record.total_time):
            record.CutCurrSig()
            tmp = record.curr_sig
            corr_list = utils_obj.CorrWith4(record,sp_list)
            toa , max_corr = utils_obj.FindTOA(corr_list,record)
            writer.writerow({"iteration":iterr,"toa_sp_1":toa[0], "toa_sp_2":toa[1], "toa_sp_3":toa[2],
                             "toa_sp_4":toa[3],"corr_1":max_corr[0],"corr_2":max_corr[1],
                             "corr_3":max_corr[2],"corr_4":max_corr[3]})
            iterr += 1

    rec_dict = utils_obj.csv2dict(record.toa_csv_path,{"iteration":[],"toa_sp_1":[],
                                                       "toa_sp_2":[], "toa_sp_3":[],
                                                       "toa_sp_4":[], "corr_1":[],"corr_2":[],
                                                       "corr_3":[],"corr_4":[]})

    print colored(" Finished parsing wav file","green")

    sp2mic = []
    # plt.figure(5)
    # plt.title('toa')
    xx = rec_dict["iteration"]
    for key in rec_dict.keys():
        if key.startswith("toa"):
            sp2mic.append(rec_dict[key])

    sp_location = utils_obj.buildSpeakersLocationMatrix(sp_list)

    if algorithm == 1:
        chan = algos.ChanAlgo()
        location_list = chan.chan_main(sp2mic, sp_location)
    elif algorithm == 2:
        taylor_obj = algos.TaylorLS()
        print "TBD"
        results_dict = {}
    elif algorithm == 3:
        LUT_obj = algos.RoomMatrix()
        location_list = LUT_obj.RoomMatMain(sp2mic,
                                            sp_list,
                                            {'x': 1.0, 'y': 1.0, 'z': 1.0},
                                            0.05,
                                            'square')
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

    if algorithm < 4:
        results_dict = utils_obj.res2csv(record.time_samples_vect,location_list,record.results_path)

    print colored('Chan_Algorithm_finished', 'green')
    for i in results_dict["Iteration"]:
        i = int(i)
        print colored("T:{0} , X:{1}, Y:{2}, Z:{3}".format(results_dict["Time [sec]"][i - 1],
                                                           results_dict["X [m]"][i - 1],
                                                           results_dict["Y [m]"][i - 1],
                                                           results_dict["Z [m]"][i - 1]), "blue")











