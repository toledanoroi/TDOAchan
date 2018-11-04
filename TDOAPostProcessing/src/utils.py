import numpy as np
import scipy.signal
import csv
import os
from scipy import signal


class UTILS(object):
    def __init__(self):


        self.yarden = []
        self.func_list = []
        self.tts_path = ''

    def _sgn(self, x):
      y = np.zeros_like(x)
      y[np.where(x >= 0)] = 1.0
      y[np.where(x < 0)] = -1.0
      return y

    def stzcr(self, x, win):
      """Compute short-time zero crossing rate."""
      if isinstance(win, str):
        win = scipy.signal.get_window(win, max(1, len(x) // 8))
      win = 0.5 * win / len(win)
      x1 = np.roll(x, 1)
      x1[0] = 0.0
      abs_diff = np.abs(self._sgn(x) - self._sgn(x1))
      return scipy.signal.convolve(abs_diff, win, mode="same")

    def ste(self, x, win):
      """Compute short-time energy."""
      if isinstance(win, str):
        win = scipy.signal.get_window(win, max(1, len(x) // 8))
      win = win / len(win)
      return scipy.signal.convolve(x ** 2, win ** 2, mode="same")

    def CorrWith4(self, speakers):
        corr_list = []
        for speaker in speakers:
            corr_list.append(np.correlate(speaker))
            corr_list.append(np.correlate(speaker.curr_sig, speaker.chirp, 'full'))
        return corr_list

    def FindTOA(self, corr_l, record, expected_signal_size, mode='max', thres=100):
        toa = []
        max_corr = []
        if mode == 'threshold':
            for v in corr_l:
                print "TBD " + str(thres)
        else:
            for v in corr_l:
                max_corr.append(max(v))
                max_index = v.argmax(axis=0)
                curr_toa = record.curr_time + (float(max_index) - expected_signal_size)/record.sample_rate
                toa.append(curr_toa)

        return toa, max_corr

    def csv2dict(self,csv_path,my_dict):

        #add if path is valid.
        with open(csv_path,'rb') as fin:
            reader = csv.DictReader(fin)
            for line in reader:
                for key,value in line.items():
                    my_dict[key].append(float(value))
        return my_dict

    def res2csv(self,time_vect,locations,path):
        with open(path, 'wb') as fout:
            writer = csv.DictWriter(fout, fieldnames=["Iteration", "Time [sec]", "X [m]", "Y [m]", "Z [m]"])
            writer.writeheader()
            iterr = 1
            for ind in range(len(time_vect)):
                writer.writerow({"Iteration": iterr, "Time [sec]":time_vect[ind],
                                 "X [m]":locations[ind][0],"Y [m]":locations[ind][1],"Z [m]":locations[ind][2]})
                iterr += 1
        results_dict = self.csv2dict(path,{"Iteration":[], "Time [sec]":[],"X [m]":[],"Y [m]":[],"Z [m]":[]})
        return results_dict

    def buildSpeakersLocationMatrix(self,sp_list):
        sp_mat = np.matrix([[sp_list[0].x, sp_list[0].y, sp_list[0].z],
                            [sp_list[1].x, sp_list[1].y, sp_list[1].z],
                            [sp_list[2].x, sp_list[2].y, sp_list[2].z],
                            [sp_list[3].x, sp_list[3].y, sp_list[3].z]])
        return sp_mat

    def CreateTTSSignal(self, msg, path='../inputs/welcome_msg.mp3'):
        from gtts import gTTS
        if os.path.isdir(path[:path.rfind('/')]):
            self.tts_path = path
        tts = gTTS(text=msg, lang='en')
        tts.save(path)
        self.Mp32Wav()
        sig_dict = self.ReadWaveFromFile()
        return sig_dict

    def Mp32Wav(self):
        from pydub import AudioSegment
        sound = AudioSegment.from_mp3(self.tts_path)
        wanted_path = self.tts_path[:-3] + 'wav'
        sound.export(wanted_path, format="wav")
        self.tts_path = wanted_path

    def ReadWaveFromFile(self, which='tts'):
        import wave
        if which == 'tts':
            tmp_path = self.tts_path
        self.sig_file_pointer = wave.open(tmp_path, 'r')
        self.num_of_channels = self.sig_file_pointer.getnchannels()
        # Extract Raw Audio from Wav File
        signal = np.fromstring(self.sig_file_pointer.readframes(-1), 'Int16')
        Fs = self.sig_file_pointer.getframerate()

        # # If Stereo
        # if self.num_of_channels == 2:
        #     print 'Our Application support only mono files'
        #     sys.exit(0)

        record_time_vector = np.linspace(0, len(signal) / Fs , num=len(signal))
        total_time = record_time_vector[-1]

        # print "Signal Parameters:\n\tFile Name:\t{0}\n\tSample rate:\t{1}" \
        #       "\n\tNumber of channels:\t{2}\n\tDuration Time:\t{3}".format(self.path,
        #                                                                    Fs,
        #                                                                    self.num_of_channels,
        #                                                                    record_time_vector[-1])

        from matplotlib import pyplot as plt
        plt.figure(1)
        plt.title(tmp_path)
        plt.plot(record_time_vector,signal)
        plt.show()
        return {'signal':signal, 'time_vect':record_time_vector, 'sample_rate':Fs, 'total_time':total_time}


class SignalHandler(object):
    def __init__(self):
        pass

    def defineParams(self,my_signal):
        self.signal = my_signal['signal']
        self.Fs = my_signal['Fs']
        self.BW = [my_signal['low_freq'], my_signal['high_freq'], abs(my_signal['high_freq'] - my_signal['low_freq'])]
        self.filtered_signal = []
        self.smoothed_signal = []
        self.new_signal = []
        self.record_time = my_signal['time_vect']

    def BPF(self, filterOrder):
        '''
        This function creates a bandpass filter with center frequency at
        centerFreq and with a passband with of pBandwidt(or pBandwidth/2 to each
        side from centerFreq). filterOrder is the order of the filter to be designed.
        Fs is the sampling frequency.
        :param filterOrder: length of filter (time domain)
        :return: update self.filtered_signal
        '''

        n = filterOrder
        pbw = (self.BW[2]) / 2
        cf = (self.BW[1] + self.BW[0]) / 2


        self.h1 = signal.firwin(numtaps=n, cutoff=[cf - pbw,cf + pbw],pass_zero=False, nyq=self.Fs / 2)
        self.mfreqz(self.h1)
        self.filtered_signal = signal.convolve(self.signal, self.h1)



    def SmoothSig(self,filterOrder,step='a'):
        '''
        smoothing the noise
        :param filterOrder:
        :param step:
        :return:
        '''
        if step == 'a':
            b, a = signal.butter(filterOrder, 0.05)
            self.smoothed_signal = signal.filtfilt(b, a, self.signal)
        elif step == 'b':
            b, a = signal.butter(filterOrder, 0.05)
            self.new_signal = signal.filtfilt(b, a, self.filtered_signal)



    def mfreqz(self, b,a=1):
        from matplotlib import pyplot as plt
        from pylab import unwrap, arctan2, imag, real, log10

        w, h = signal.freqz(b,a)
        h_dB = 20 * log10(abs(h))
        plt.subplot(211)
        plt.plot(w/max(w), h_dB)
        plt.ylim(-150, 5)
        plt.ylabel('Magnitude (db)')
        plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
        plt.title(r'Frequency response')
        plt.subplot(212)
        h_Phase = unwrap(arctan2(imag(h),real(h)))
        plt.plot(w/max(w),h_Phase)
        plt.ylabel('Phase (radians)')
        plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
        plt.title(r'Phase response')
        plt.subplots_adjust(hspace=0.5)
        plt.show()

    def Plotsignaldiff(self,mode='f'):
        old = self.signal

        if mode == 'f':
            new = self.filtered_signal
        elif mode == 's':
            new = self.smoothed_signal
        elif mode == 'a':
            new1 = self.filtered_signal
            new2 = self.smoothed_signal
            new3 = self.new_signal

        from matplotlib import pyplot as plt

        fig = plt.figure()
        plt.plot(self.record_time,old,)



if __name__ == '__main__':
    a = UTILS()
    sig_d = a.CreateTTSSignal('Welcome to the ultrasonic sound localization system by Roi Toledano and Yarden Avraham.')
    for item in sig_d.items:
        print item

    from src.wave2toa import Speaker
    from src.wave2toa import recwav
    rec = recwav()
    rec.change_path('../inputs/second_test.wav','in')
    rec.PlotSignal('second.wav')
    s = Speaker()
    s.Define_ID(1)
    s.BuildChirp()
    sp_sig = {'signal': s.chirp, 'Fs': 88200, 'low_freq': 27500, 'high_freq': 30000, 'mat_signal': s.matlab_chirp}
    rec_sig = {'signal': rec.signal, 'Fs': 88200, 'low_freq': 7500, 'high_freq': 30000, 'time': rec.rec_time}
    b = SignalHandler()
    c = SignalHandler()
    b.defineParams(sp_sig)
    c.defineParams(rec_sig)

    b.BPF(139)
    c.BPF(139)


    # b.SmoothSig(139, step='a')
    # b.SmoothSig(139, step='b')
    # c.SmoothSig(139, step='a')
    # c.SmoothSig(139, step='b')

    print "Finish All"