import numpy as np
import scipy.signal
import csv
import os
from scipy import signal
import statistics as stats
from math import exp
from matplotlib import pyplot as plt


class UTILS(object):
    def __init__(self):
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

    def CreateTDOAlistBySPx(self,sp_number,data_dict):
        TDOA_for_outliers = {"tdoa_sp_1": [data_dict['toa_sp_1'][i] - data_dict['toa_sp_' + str(sp_number)][i] for i in range(len(data_dict['toa_sp_1']))],
                             "tdoa_sp_2": [data_dict['toa_sp_2'][i] - data_dict['toa_sp_' + str(sp_number)][i] for i in range(len(data_dict['toa_sp_2']))],
                             "tdoa_sp_3": [data_dict['toa_sp_3'][i] - data_dict['toa_sp_' + str(sp_number)][i] for i in range(len(data_dict['toa_sp_3']))],
                             "tdoa_sp_4": [data_dict['toa_sp_4'][i] - data_dict['toa_sp_' + str(sp_number)][i] for i in range(len(data_dict['toa_sp_4']))]
                             }
        return TDOA_for_outliers

    def Sp2MicToTDOA(self,sp2mic):
        from collections import OrderedDict as OD
        rows = len(sp2mic)
        measuredTDOA_vectors = OD()
        for i in range(rows):
            tmpvector = np.transpose(np.array(sp2mic) - np.array(sp2mic[i]))
            measuredTDOA_vectors['TDOA_from_sp' + str(i+1)] = tmpvector
        return measuredTDOA_vectors

    def AveragingSamples(self,avg_list,time_vect,avg_dim,TDOA_vect):
        from statistics import mean
        from collections import OrderedDict as OD
        v = [[],[],[],[]]
        time_avg = []
        v_avg = OD()
        v_avg['TDOA_from_sp1'] = []
        v_avg['TDOA_from_sp2'] = []
        v_avg['TDOA_from_sp3'] = []
        v_avg['TDOA_from_sp4'] = []

        last = 0
        last_ = 0
        curr_ = 0
        curr = 0
        for avg in avg_list:
            curr += avg
            curr_ += avg_dim
            for i in range(len(v)):
                v[i] = TDOA_vect['TDOA_from_sp' + str(i+1)][last:curr]
                v_avg['TDOA_from_sp' + str(i+1)].append([mean([k[j] for k in v[i]]) for j in range(len(v[i][0]))])
            time_avg.append(mean(time_vect[last_:curr_]))
            last = curr
            last_ = curr_

        return v_avg, time_avg



    def CorrSpeakerSig(self, speaker):
        a = np.correlate(speaker.proccessed_signal.filtered_signal, speaker.matlab_chirp, 'full')
        b = np.correlate(speaker.proccessed_signal.signal, speaker.matlab_chirp, 'full')
        return a, b

    def CorrWith4(self, speakers):
        corr_list = []
        for speaker in speakers:
            corr_list.append(np.correlate(speaker.proccessed_signal.filtered_signal, speaker.chirp, 'full'))
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

    def res2csv(self,locations,path):
        with open(path, 'wb') as fout:
            writer = csv.DictWriter(fout, fieldnames=["Iteration", "Time [sec]", "X [m]", "Y [m]", "Z [m]", "Error"])
            writer.writeheader()
            iterr = 1
            for locat in locations:
                writer.writerow({"Iteration": iterr, "Time [sec]":locat[0],
                                 "X [m]":locat[1][0],"Y [m]":locat[1][1],"Z [m]":locat[1][2],"Error":locat[2]})
                iterr += 1
        results_dict = self.csv2dict(path,{"Iteration": [], "Time [sec]": [], "X [m]": [], "Y [m]": [], "Z [m]": [], "Error": []})
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

    def CalcGaussianParams_nD(self,vectors):
        mean_list = []
        std_list = []
        variance_list = []
        median_list = []
        for vec in vectors:
            mean_list.append(stats.mean(vec))
            std_list.append(stats.stdev(vec))
            variance_list.append(stats.variance(vec))
            median_list.append(stats.median(vec))

        return {'mean': mean_list,
                'stdev': std_list,
                'variance': variance_list,
                'median': mean_list
                }

    def GaussianFunc(self,sample,statistics_dict):
        if len(sample) != len(statistics_dict['mean']):
            raise ValueError()
        frac = 0
        for i in range(len(sample)):
            frac += self.ExpFrac(sample[i],statistics_dict['mean'][i],statistics_dict['variance'][i])
        return exp(-1 * frac)

    def ExpFrac(self,x,mu,var):
        return ((x-mu) ** 2) / (2 * var)

    def Find_Outliers(self,vectors):
        gauss_prm = self.CalcGaussianParams_nD(vectors)
        is_outlier = []
        for i in range(len(vectors[0])):
            f = self.GaussianFunc([v[i] for v in vectors],gauss_prm)
            if f < exp(-1.5):                                   # distance of 3 sigmas in all dimentions
                is_outlier.append(True)
            else:
                is_outlier.append(False)

        return is_outlier


    def Throw_Outliers(self,dict_r, dict_im, dim, keys=["toa_sp_1","toa_sp_2","toa_sp_3","toa_sp_4"]):
        '''
        This function cut samples to group in dim size , search for outliers and throw them from the data set.

        :param dict_r: TOA samples dictionary
        :param dict_im: TDOA samples dictionary
        :param dim: samples groups size
        :param keys: key for the new dictionary
        :return: new TOA samples dictionary and list of the new groups sizes to average in the algorithm
        '''
        rec_dict1 = {keys[0]: [],
                     keys[1]: [],
                     keys[2]: [],
                     keys[3]: []
                     }
        to_average = []
        v_r = [[],[],[],[]]
        v = [[],[],[],[]]
        for i in range(len(dict_im['tdoa_sp_1'])):
            if len(v[3]) == 5 & len(v[1]) == 5 & len(v[2]) == dim:
                outliers_list = self.Find_Outliers(v[1:])
                c = 0
                for k in range(len(outliers_list)):
                    if outliers_list[k] == False:
                        c += 1
                        for j in range(len(v_r)):
                            rec_dict1['toa_sp_' + str(j+1)].append(v_r[j][k])
                to_average.append(c)
                v_r = [[], [], [], []]
                v = [[], [], [], []]

            for l in range(len(v)):
                v[l].append(dict_im['tdoa_sp_' + str(l + 1)][i])
            for h in range(len(v_r)):
                v_r[h].append(dict_r['toa_sp_' + str(h+1)][i])

            # add exception for the reminder


        return rec_dict1 , to_average

    def ScatterPlot3D(self,x_pnt,y_pnt,z_pnt,title,labels,limits):
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # ax = Axes3D(fig)
        # For each set of style and range settings, plot n random points in the box
        # # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
        # for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
        for i in range(len(x_pnt)):
            ax.scatter(x_pnt[i], y_pnt[i], z_pnt[i], c='b', marker='o')

        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_zlabel(labels[2])
        ax.set_title(title)
        Axes3D.set_xlim(ax, right=limits[0][1], left=limits[0][0])
        Axes3D.set_ylim(ax, bottom=limits[1][0], top=limits[1][1])
        Axes3D.set_zlim(ax, bottom=limits[2][0], top=limits[2][1])

        plt.show()

    def ClusteringnD(self,dim,data,time_vect):
        '''
        The function get data sets with dimention dim, generate vectors and search for groups with same parameters
        :param dim: The dimention of one vector
        :param data: vectors set
        :return: dicitionary of groups groups = { 'g_1' : [<sample number> ,<timestamp>, <vector>], ... }
        '''
        pass


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

    def BPF(self, filterOrder, plotting=False):
        '''
        This function creates a bandpass filter with center frequency at
        centerFreq and with a passband with of pBandwidt(or pBandwidth/2 to each
        side from centerFreq). filterOrder is the order of the filter to be designed.
        Fs is the sampling frequency.
        :param filterOrder: length of filter (time domain)
        :return: update self.filtered_signal
        '''

        n = filterOrder
        pbw = 1.1 * ((self.BW[2]) / 2)
        cf = (self.BW[1] + self.BW[0]) / 2


        self.h1 = signal.firwin(numtaps=n, cutoff=[cf - pbw,cf + pbw],pass_zero=False, nyq=self.Fs / 2)
        # if plotting:
        #     self.mfreqz(self.h1)
        self.filtered_signal = signal.convolve(self.signal, self.h1)
        self.record_time_with_filter_delay = np.linspace(0, float(len(self.filtered_signal)) / self.Fs, num=len(self.filtered_signal))



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
        plt.grid()
        plt.ylim(-150, 5)
        plt.ylabel('Magnitude (db)')
        plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
        plt.title(r'Frequency response')
        plt.subplot(212)
        h_Phase = unwrap(arctan2(imag(h),real(h)))
        plt.plot(w/max(w),h_Phase)
        plt.grid()
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
    rec.change_path('../inputs/blackmanharris5ms/1.WAV','in')
    rec.PlotSignal('blackmanharris5ms.wav')
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