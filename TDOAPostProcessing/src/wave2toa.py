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



def resampling(speaker,rec_Fs):
    factor = int((rec_Fs/250000.0)*len(speaker.matlab_chirp))
    f = signal.resample(speaker.matlab_chirp,factor)
    x = np.linspace(0, 0.001, 251, endpoint=False)
    x1 = np.linspace(0, 0.001, factor, endpoint=False)
    plt.plot(x, speaker.matlab_chirp, 'go-', x1, f, '.-',x1,speaker.chirp,'.-', 1, speaker.matlab_chirp[0], 'ro')
    plt.legend(['data', 'resampled','python generated'], loc='best')
    plt.show()

    fft1   = fftpkt.fftshift(fftpkt.fft(speaker.matlab_chirp))
    freqs  = fftpkt.fftshift(fftpkt.fftfreq(len(speaker.matlab_chirp), 1.0 / 88200))
    fft22  = fftpkt.fftshift(fftpkt.fft(f))
    freqs2 = fftpkt.fftshift(fftpkt.fftfreq(len(f), 1.0 / 88200))
    fft33  = fftpkt.fftshift(fftpkt.fft(speaker.chirp))
    freqs3 = fftpkt.fftshift(fftpkt.fftfreq(len(speaker.chirp), 1.0 / 88200))

    plt.figure(10)
    plt.plot(freqs, abs(fft1), 'r')
    plt.plot(freqs2, abs(fft22), 'b')
    plt.plot(freqs3, abs(fft33), 'g')
    plt.legend(['matlab 250000', 'resampled', 'python generated'], loc='best')
    plt.show()



class recwav():
    def __init__(self):
        self.path = '../inputs/first_test.wav'
        self.sample_rate = 44100
        self.toa_csv_path = '../output/toa_record.csv'
        self.record_sig = []
        self.num_of_channels = 1


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

        print "Signal Parameters:\n\tFile Name:\t{0}\n\tSample rate:\t{1}" \
              "\n\tNumber of channels:\t{2}\n\tDuration Time:\t{3}".format(self.path,
                                                                           self.sample_rate,
                                                                           self.num_of_channels,
                                                                           self.rec_time[-1])

        plt.figure(1)
        plt.title('Signal Wave: ' + title )
        plt.plot(self.rec_time,self.signal)
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

    def change_path(self,path,mode = 'in'):
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

        # plt.figure(3)
        # plt.title('one cycle')
        # plt.plot(tmprec_time,tmpsig)
        # plt.show()

        plt.figure(4)
        plt.title('one cycle FFT:')
        plt.plot(tmpfreqs,abs(tmpfft))
        plt.grid()
        plt.show()

    def FindRXTimestamps(self):
        return 0

    def Spectogram(self):
        f,t,Sxx = signal.spectrogram(self.signal,self.sample_rate, nperseg=64)
        plt.pcolormesh(t,f,(10 * np.log10(Sxx)), cmap='inferno')
        plt.colorbar()
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        plt.show()


class Speaker():
    def __init__(self):
        self.id = 555
        self.x = 0
        self.y = 0
        self.z = 0
        self.chirp = []

    def BuildChirp(self):
        Fs = 88200
        tt = np.linspace(0,0.001,Fs*0.001)
        matlabchirps = io.loadmat('../inputs/all_chirp.mat')
        chirps = matlabchirps['allchirp']
        if self.id == 1:
            tmpchirp = signal.chirp(tt,30000,0.001,27500)
            nrm = np.linspace(1,0.02,num=len(tmpchirp))
            return np.multiply(nrm,tmpchirp) , chirps[:,0]
        if self.id == 2:
            tmpchirp = signal.chirp(tt,15000,0.001,12500)
            nrm = np.linspace(1,0.02,num=len(tmpchirp))
            return np.multiply(nrm,tmpchirp) , chirps[:,1]
        if self.id == 3:
            tmpchirp = signal.chirp(tt,20000,0.001,17500)
            nrm = np.linspace(1,0.02,num=len(tmpchirp))
            return np.multiply(nrm,tmpchirp) , chirps[:,2]
        if self.id == 4:
            tmpchirp = signal.chirp(tt,25000,0.001,22500)
            nrm = np.linspace(1,0.02,num=len(tmpchirp))
            return np.multiply(nrm,tmpchirp) , chirps[:,3]

    def Define_ID(self,my_id):
        self.id = my_id
        self.chirp, self.matlab_chirp = self.BuildChirp()


    def DefineLocation(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z




if __name__ == '__main__':
    sp1 = Speaker()
    sp1.Define_ID(1)
    sp2 = Speaker()
    sp2.Define_ID(2)
    sp3 = Speaker()
    sp3.Define_ID(3)
    sp4 = Speaker()
    sp4.Define_ID(4)


    record = recwav()
    record.change_path('../inputs/second_test.wav','in')
    record.PlotSignal('second.wav')
    #record.PlotFFT('second.wav')
    # record.plotOnlyOneCycle()
    record.Spectogram()

    corr1 = np.correlate(record.signal,sp1.matlab_chirp, 'full')
    corr2 = np.correlate(record.signal, sp2.matlab_chirp, 'full')
    corr3 = np.correlate(record.signal, sp3.matlab_chirp, 'full')
    corr4 = np.correlate(record.signal, sp4.matlab_chirp, 'full')


    resampling(sp1,record.sample_rate)


    time1 = np.linspace(0, len(corr1) / record.sample_rate, num=len(corr1))
    plt.figure(5)
    plt.title('correlation')
    plt.plot(time1, abs(corr1),'r')
    plt.plot(time1, abs(corr2), 'b')
    # plt.plot(time1, abs(corr3), 'g')
    # plt.plot(time1, abs(corr4), 'p')
    plt.grid()
    plt.show()


    # f, (ax1, ax2,ax3,ax4) = plt.subplots(2, 2, sharey=True)
    # ax1.plot(record.rec_time, corr1)
    # ax1.set_title('chirp 1')
    # ax2.plot(record.rec_time, corr2)
    # ax2.set_title('chirp 2')
    # ax3.plot(record.rec_time, corr3)
    # ax3.set_title('chirp 3')
    # ax4.plot(record.rec_time, corr4)
    # ax4.set_title('chirp 4')

    # fig = plt.figure()
    #
    # plt.subplot(2, 2, 1)
    # plt.plot(record.rec_time, corr1)
    # plt.title('1')
    #
    # plt.subplot(2, 2, 2)
    # plt.plot(record.rec_time, corr2)
    # plt.title('2')
    #
    # plt.subplot(2, 2, 3)
    # plt.plot(record.rec_time, corr3)
    # plt.title('3')
    #
    # plt.subplot(2, 2, 4)
    # plt.plot(record.rec_time, corr4)
    # plt.title('4')
    #
    # plt.show()








