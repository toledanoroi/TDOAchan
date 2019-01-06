# project: TDOA based ultrasonic sound localization system
# authors: Roi Toledano & Yarden Avraham
# lab    : Bats lab
# guide  : PhD Yossi Yovel

import numpy as np
import pandas as pd
import scipy.signal
import csv
import os
from scipy import signal
import statistics as stats
from math import exp
from matplotlib import pyplot as plt
from numpy import concatenate
from scipy.spatial import ConvexHull


class UTILS(object):
    def __init__(self):
        self.func_list = []
        self.tts_path = ''

    def _sgn(self, x):
        '''
        returns the sign of each number in the signal
        :param x: signal
        :return: list of 1/-1 the describes the sign for each index in the signal
        '''
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

    def CreateTDOAlistBySPx(self, sp_number, data_dict):
        '''
        calculate TDOA matrix from toa measurements , TDOA from chosen sp_number
        :param sp_number: number of speaker to calculate TDOA from him.
        :param data_dict: dictionary of all the toa measurements.
        :return: TDOA_for_outliers , dictionary of TDOA array
        '''
        TDOA_for_outliers = {"tdoa_sp_1": [data_dict['toa_sp_1'][i] - data_dict['toa_sp_' + str(sp_number)][i] for i in range(len(data_dict['toa_sp_1']))],
                             "tdoa_sp_2": [data_dict['toa_sp_2'][i] - data_dict['toa_sp_' + str(sp_number)][i] for i in range(len(data_dict['toa_sp_2']))],
                             "tdoa_sp_3": [data_dict['toa_sp_3'][i] - data_dict['toa_sp_' + str(sp_number)][i] for i in range(len(data_dict['toa_sp_3']))],
                             "tdoa_sp_4": [data_dict['toa_sp_4'][i] - data_dict['toa_sp_' + str(sp_number)][i] for i in range(len(data_dict['toa_sp_4']))]
                             }
        return TDOA_for_outliers

    def Sp2MicToTDOA(self, sp2mic):
        '''
        calculate TDOA matrix from sp2mic toa measurements in regard of all the speakers
        :param sp2mic: list of all the toa measurements.
        :return:
        '''
        from collections import OrderedDict as OD
        rows = len(sp2mic)
        measuredTDOA_vectors = OD()
        for i in range(rows):
            tmpvector = np.transpose(np.array(sp2mic) - np.array(sp2mic[i]))
            measuredTDOA_vectors['TDOA_from_sp' + str(i+1)] = tmpvector
        return measuredTDOA_vectors

    def AveragingSamples(self, avg_list, time_vect, avg_dim, TDOA_vect):
        '''
        calculate average of tdoa samples of all average groups (after outliers throwing process)
        :param avg_list: list of amounts of samples in each group to average.
        :param time_vect: samples timestamps
        :param avg_dim: The group size before outliers throwing process
        :param TDOA_vect: The TDOA samples
        :return:
        '''
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
        '''
        correlation between filtered signal and expected matlab generated signal
        and correlation between unfiltered signal and expected signal
        :param speaker: the speaker to correlate
        :return: 2 correlation vectors
        corr_after_filter, corr_before_filter
        '''
        corr_after_filter = np.correlate(speaker.proccessed_signal.filtered_signal, speaker.matlab_chirp, 'full')
        corr_before_filter = np.correlate(speaker.proccessed_signal.signal, speaker.matlab_chirp, 'full')
        return corr_after_filter, corr_before_filter

    def csv2dict(self, csv_path, my_dict={}):
        '''
        reads csv and generate/append to dictionary
        :param csv_path: the csv path
        :param my_dict: the dictionary to append ,
        default: generation my_dict = {}
        :return:
        '''

        #add if path is valid.
        if os.path.isfile(csv_path):
            with open(csv_path,'rb') as fin:
                reader = csv.DictReader(fin)
                for line in reader:
                    for key,value in line.items():
                        my_dict[key].append(float(value))
        return my_dict

    def res2csv(self, locations, path, error, ex, ey, ez, pnt_name):
        '''
        takes results dictionary and generate csv of it.
        :param locations: results dictionary
        :param path: path of results csv
        :return:
        '''
        with open(path, 'wb') as fout:
            writer = csv.DictWriter(fout, fieldnames=["Iteration", "Time [sec]", "X [m]", "Y [m]", "Z [m]", "cost_function_value","point_set",'E[x]','E[y]','E[z]','Error'])
            writer.writeheader()
            iterr = 1
            for locat in locations:
                writer.writerow({"Iteration": iterr, "Time [sec]":locat[0],
                                 "X [m]": locat[1][0],
                                 "Y [m]": locat[1][1],
                                 "Z [m]": locat[1][2],
                                 "cost_function_value": locat[2],
                                 "point_set": pnt_name,
                                 'E[x]': ex,
                                 'E[y]': ey,
                                 'E[z]': ez,
                                 'Error': error
                                 }
                                )
                iterr += 1
        results_dict = self.csv2dict(path,{"Iteration": [],
                                           "Time [sec]": [],
                                           "X [m]": [],
                                           "Y [m]": [],
                                            "Z [m]": [],
                                            "cost_function_value": [],
                                            "point_set": [],
                                            'E[x]': [],
                                            'E[y]': [],
                                            'E[z]': [],
                                            'Error': []
                                           })
        return results_dict

    def buildSpeakersLocationMatrix(self,sp_list):
        '''
        generate locations matrix
        |-------------------|
        |  x1     y1     z1 |
        |  x2     y2     z2 |
        |  x3     y3     z3 |
        |  x4     y4     z4 |
        |-------------------|
        :param sp_list: speakers objects list
        :return:
        '''
        sp_mat = np.matrix([[sp_list[0].x, sp_list[0].y, sp_list[0].z],
                            [sp_list[1].x, sp_list[1].y, sp_list[1].z],
                            [sp_list[2].x, sp_list[2].y, sp_list[2].z],
                            [sp_list[3].x, sp_list[3].y, sp_list[3].z]])
        return sp_mat

    def CreateTTSSignal(self, msg, path='../inputs/welcome_msg.mp3'):
        '''
        generate text to speech message mp3 file, convert it to wav file and returns it as ndarray
        :param msg: text to speech message
        :param path: path for new mp3 file.
        default: path='../inputs/welcome_msg.mp3'
        :return: signal taps and parameters dictionary
        {'signal':signal, 'time_vect':record_time_vector, 'sample_rate':Fs, 'total_time':total_time}
        '''
        from gtts import gTTS
        if os.path.isdir(path[:path.rfind('/')]):
            self.tts_path = path
        tts = gTTS(text=msg, lang='en')
        tts.save(path)
        self.Mp32Wav()
        sig_dict = self.ReadWaveFromFile()
        return sig_dict

    def Mp32Wav(self):
        '''
        converts mp3 file to wav file.
        :return: path for new .wav file
        '''
        from pydub import AudioSegment
        sound = AudioSegment.from_mp3(self.tts_path)
        wanted_path = self.tts_path[:-3] + 'wav'
        sound.export(wanted_path, format="wav")
        self.tts_path = wanted_path

    def ReadWaveFromFile(self, which='tts'):
        '''
        read wav file and converts it to ndarray.
        finds all parameters of the signal and returns it as dictionary
        :param which: type of signal
        :return: dictionary of signal parameters
        '''
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
        plt.show(block=True)
        return {'signal':signal, 'time_vect':record_time_vector, 'sample_rate':Fs, 'total_time':total_time}

    def CalcGaussianParams_nD(self, vectors):
        '''
        calculate gaussian parameters from samples
        mean , std , variance
        :param vectors: samples vectors to calculate the gaussian parameters
        :return: dictionary with gaussian parameters
        {'mean': mean_list,
                'stdev': std_list,
                'variance': variance_list,
                'median': mean_list }
        '''
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

    def GaussianFunc(self, sample, statistics_dict):
        '''
        Calculate Gaussian function value of sample
        :param sample:
        :param statistics_dict: gaussian parameters
        :return: Gaussian function value
        '''
        if len(sample) != len(statistics_dict['mean']):
            raise ValueError()
        frac = 0
        for i in range(len(sample)):
            frac += self.ExpFrac(sample[i],
                                 statistics_dict['mean'][i],
                                 statistics_dict['variance'][i])
        return exp(-1 * frac)

    def ExpFrac(self, x, mu, var):
        '''
        Calculate the Gaussian exponent argument for one dimension
        :param x: the sample
        :param mu: E[x]
        :param var: Variance[x]
        :return: Gaussian exponent argument
        '''
        if (var == 0):
            var = 0.001
        return ((x-mu) ** 2) / (2 * var)

    def Find_Outliers(self, vectors):
        '''
        search for outliers in group of vectors according to Gaussian nD model
        :param vectors: group of samples
        :return: Boolean list of which sample is outliers and which not.
        True -> not outlier , False-> outlier
        '''
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

    def ScatterPlot3D(self, x_pnt, y_pnt, z_pnt, title, labels, limits, cvx1=None, cvx2=None, expected=None):
        '''
        plotting 3D scatter plot of locations in the room , plot the room edges too.
        :param x_pnt: list of x samples
        :param y_pnt: list of y samples
        :param z_pnt: list of z samples
        :param title: The plot title
        :param labels: axis labels [<first>, <second>, <third>]
        :param limits: the max size in the plot for each axis
        :param cvx1: the convex hull of the room
        default = None
        :param cvx2: the convex hull of unallowed volume in the room
        default = None
        :param expected: expected location
        default = None
        '''
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
        if cvx1 is not None:
            pts = cvx1.points
            ax.plot(pts.T[0], pts.T[1], pts.T[2], "ko")
            for s in cvx1.simplices:
                s = np.append(s, s[0])  # Here we cycle back to the first coordinate
                ax.plot(pts[s, 0], pts[s, 1], pts[s, 2], "r-")

        if cvx2 is not None:
            pts2 = cvx2.points
            ax.plot(pts2.T[0], pts2.T[1], pts2.T[2], "mo")
            for s in cvx2.simplices:
                s = np.append(s, s[0])  # Here we cycle back to the first coordinate
                ax.plot(pts2[s, 0], pts2[s, 1], pts2[s, 2], "m-")

        if expected is not None:
            for pnt in expected:
                ax.scatter(pnt[0], pnt[1], pnt[2], c='g', marker='*')

        plt.show(block=True)

    def PlotCvxHull(self, hull, cvx2=None):
        '''
        Plot convex hull --> for debug
        :param hull: the convex hull of the room
        :param cvx2: the convex hull of unallowed volume in the room
        default = None
        '''
        from mpl_toolkits.mplot3d import Axes3D
        from scipy.spatial import ConvexHull

        # 8 points defining the cube corners
        pts = hull.points
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

        # Plot defining corner points
        ax.plot(pts.T[0], pts.T[1], pts.T[2], "ko")

        # 12 = 2 * 6 faces are the simplices (2 simplices per square face)
        for s in hull.simplices:
            s = np.append(s, s[0])  # Here we cycle back to the first coordinate
            ax.plot(pts[s, 0], pts[s, 1], pts[s, 2], "r-")

        if cvx2 is not None:
            pts2 = cvx2.points
            ax.plot(pts2.T[0], pts2.T[1], pts2.T[2], "mo")
            for s_ in cvx2.simplices:
                s_ = np.append(s_, s_[0])  # Here we cycle back to the first coordinate
                ax.plot(pts2[s_, 0], pts2[s_, 1], pts2[s_, 2], "m-")
        # # Make axis label
        # for i in ["x", "y", "z"]:
        #     ["myconvexhull"][1]
        #     eval("ax.set_{:s}label('{:s}')".format(i, i))

        plt.show(block=True)

    def ScatterPlot2D(self, x_pnt, y_pnt, title, labels, limits, cvx1=None, cvx2=None, expected=None):
        '''
        plotting 2D scatter plot of locations in the room , plot the room edges too.
        :param x_pnt: list of x samples
        :param y_pnt: list of y samples
        :param title: The plot title
        :param labels: axis labels [<first>, <second>, <third>]
        :param limits: the max size in the plot for each axis
        :param cvx1: the convex hull of the room
        default = None
        :param cvx2: the convex hull of unallowed volume in the room
        default = None
        :param expected: expected location
        default = None
        '''
        fig = plt.figure()
        c = 0
        ax = fig.add_subplot(111)
        for i in range(len(x_pnt)):
            ax.scatter(x_pnt[i], y_pnt[i], c='b', marker='o')
            if np.mod(i, 3) == 0:
                c += 1
            #     if (c == 3):
            #         c += 1
            #     elif (c == 12):
            #         c += 1
            # ax.annotate(str(c), (x_pnt[i], y_pnt[i]))


        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_title(title)
        ax.set_xlim(right=limits[0][1], left=limits[0][0])
        ax.set_ylim(bottom=limits[1][0], top=limits[1][1])

        if cvx1 is not None:
            pts = cvx1.points
            ax.plot(pts.T[0], pts.T[1], "ko")
            for s in cvx1.simplices:
                s = np.append(s, s[0])
                ax.plot(pts[s, 0], pts[s, 1], "r-")

        if cvx2 is not None:
            pts2 = cvx2.points
            ax.plot(pts2.T[0], pts2.T[1], "mo")
            for s in cvx2.simplices:
                s = np.append(s, s[0])  # Here we cycle back to the first coordinate
                ax.plot(pts2[s, 0], pts2[s, 1], "m-")

        if expected is not None:
            for ii, pnt in enumerate(expected):
                ax.scatter(pnt[0], pnt[1], c='g', marker='*')
                # if np.mod(ii, 3) == 0:
                #     ax.annotate('a' + str(1 + ii/2), (pnt[0], pnt[1]))


        plt.xticks(np.arange(limits[0][0], limits[0][1], step=0.3))
        plt.yticks(np.arange(limits[1][0], limits[1][1], step=0.3))
        plt.grid()
        # plt.legend(['room edge points', 'room walls', 'not room points', 'non room walls egdes', 'algorithm results', 'expected'])
        plt.show()

    def CalcErrorFromExpected(self, pnt, expected):
        '''
        Calculate the euclidean distance between 2 ndim points
        :param pnt: the tested point
        :param expected: the expected results point (real location 3D/2D)
        :return: distance
        :type float64
        '''
        _pnt = np.array(pnt)
        _expected = np.array(expected)
        differ = (abs(_pnt - _expected)) ** 2
        return np.sqrt(differ.sum(-1))

    def ClusteringnD(self,dim,data,time_vect):
        '''
        The function get data sets with dimention dim, generate vectors and search for groups with same parameters
        :param dim: The dimention of one vector
        :param data: vectors set
        :return: dicitionary of groups groups = { 'g_1' : [<sample number> ,<timestamp>, <vector>], ... }
        '''
        from sklearn.cluster import MeanShift, estimate_bandwidth
        # from sklearn.datasets.samples_generator import make_blobs
        bandwidth = estimate_bandwidth(data, quantile=0.2, n_samples=500)

        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(data)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_

        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)

        print("number of estimated clusters : %d" % n_clusters_)
        pass

    def DefineRoom(self, room_list, not_edge_dots, shapewithnonedgepoint=True, plotting=False):
        '''
        generate convex hulls of the room from points 3D and 2D
        :param room_list: 3D arrays array of all edges in the room
        :param not_edge_dots: 3D arrays array of all the not allowed volume in the room
        :param shapewithnonedgepoint: boolean , defines if there are not allowed volumes in the room
        :param plotting: boolean, defines if we want to plot the convex hull or not.
        :return:
        '''
        from numpy import array
        from scipy.spatial import ConvexHull
        room3D = array(room_list)
        list_2d = [array([pnt[0], pnt[1]]) for pnt in room_list]
        room2D = array(list_2d)
        if not shapewithnonedgepoint:

            triangle3D = array(not_edge_dots)
            list_2d_non = [array([pnt[0], pnt[1]]) for pnt in not_edge_dots]
            triangle2D = array(list_2d_non)

            non_hull = ConvexHull(triangle3D)
            non_hull2d = ConvexHull(triangle2D)

            hull = ConvexHull(room3D)
            hull2d = ConvexHull(room2D)

            if plotting:
             self.PlotCvxHull(hull,  non_hull)

            return hull, non_hull, hull2d, non_hull2d

        hull = ConvexHull(room3D)
        hull2d = ConvexHull(room2D)

        return hull, None, hull2d, None

    def IsPntInConvexHull(self,hull, pnt):
        '''
        Checks if `pnt` is inside the convex hull.
        :param `hull` -- a QHull ConvexHull object
        :param `pnt` -- point array of shape (3,)
        '''
        new_hull = ConvexHull(concatenate((hull.points, [pnt])))
        if np.array_equal(new_hull.vertices, hull.vertices):
            return True
        return False

    def MergecsvAndGenerateForPlotting(self,folder,prefix='a'):
        import fnmatch
        import os
        err2d = []
        err3d = []
        matches = []
        # a = os.walk('/Users/roitoledano/Desktop/WLAN_vids')
        # for dirs in a:
        #     print dirs
        for root, dirnames, filenames in os.walk(folder):
            for filename in fnmatch.filter(filenames, prefix + '*.csv'):
                match = os.path.join(root, filename)
                matches.append(pd.read_csv(match))




class SignalHandler(object):
    def __init__(self):
        pass

    def defineParams(self, my_signal):
        '''
        define all SignalHandler object parameters
        :param my_signal: signal parameters
        '''
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
        #
        #     self.mfreqz(self.h1)
        self.filtered_signal = signal.convolve(self.signal, self.h1)
        self.record_time_with_filter_delay = np.linspace(0, float(len(self.filtered_signal)) / self.Fs, num=len(self.filtered_signal))

    def SmoothSig(self,filterOrder,step='a'):
        '''
        smoothing the noise
        we don't use it in the code, maybe needed
        :param filterOrder:
        :param step:
        '''
        if step == 'a':
            b, a = signal.butter(filterOrder, 0.05)
            self.smoothed_signal = signal.filtfilt(b, a, self.signal)
        elif step == 'b':
            b, a = signal.butter(filterOrder, 0.05)
            self.new_signal = signal.filtfilt(b, a, self.filtered_signal)

    def mfreqz(self, b,a=1):
        '''
        plotting freqz of filter , like matlab representation.
        :param b: nominator
        :param a: denominator
        default: a = 1
        '''
        from matplotlib import pyplot as plt
        from pylab import unwrap, arctan2, imag, real, log10

        w, h = signal.freqz(b,a)
        h_dB = 20 * log10(abs(h))
        plt.subplot(211)
        plt.plot(w/max(w), h_dB)
        plt.grid()
        plt.ylim(-150, 5)
        plt.ylabel('Magnitude (dB)')
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
        plt.show(block=True)

        plt.figure()
        plt.plot(b, 'bo-', linewidth=2)
        plt.title('Filter Coefficients (%d taps)' % len(b))
        plt.grid(True)
        plt.show()

    def Plotsignaldiff(self,mode='f'):
        '''
        plot the signal to show the difference between operation
        :param mode: {'f' : only filtered_signal , 's': only smoothed_signal, 'a': all signals}
        :return:
        '''
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
        plt.plot(self.record_time, old)
        if mode == 'a':
            plt.plot(self.record_time, new1[:len(old)])
            plt.plot(self.record_time, new2[:len(old)])
            plt.plot(self.record_time, new3[:len(old)])
        else:
            plt.plot(self.record_time, new[:len(old)])
        plt.show()


# only for debug
# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    a = UTILS()
    # sig_d = a.CreateTTSSignal('Welcome to the ultrasonic sound localization system by Roi Toledano and Yarden Avraham.')
    # for item in sig_d.items:
    #     print item


    room3D = np.array([np.array([0, 0, 0]),
                       np.array([0, 0, 2.4]),
                       np.array([3.9, 0, 0]),
                       np.array([3.9, 0, 2.4]),
                       np.array([3.9, 4.04, 0]),
                       np.array([3.9, 4.04, 2.4]),
                       np.array([0.6, 4.04, 0]),
                       np.array([0.6, 4.04, 2.4]),
                       np.array([0.6, 1.04, 0]),
                       np.array([0.6, 1.04, 2.4]),
                       np.array([0, 1.04, 0]),
                       np.array([0, 1.04, 2.4]),
                       ])
    triangle3D = np.array([np.array([0.6, 4.04, 0]),
                           np.array([0.6, 4.04, 2.4]),
                           np.array([0.6, 1.04, 0]),
                           np.array([0.6, 1.04, 2.4]),
                           np.array([0, 1.04, 0]),
                           np.array([0, 1.04, 2.4]),
                           ])


    room2D = np.array([np.array([0, 0]),
                       np.array([3.9, 0]),
                       np.array([3.9, 4.04]),
                       np.array([0.6, 4.04]),
                       np.array([0.6, 1.04]),
                       np.array([0, 1.04])
                       ])
    triangle2D = np.array([np.array([0.6, 4.04]),
                           np.array([0.6, 1.04]),
                           np.array([0, 1.04])
                           ])


    hull2d = ConvexHull(room2D)
    non_hull2d = ConvexHull(triangle2D)
    hull = ConvexHull(room3D)
    non_hull = ConvexHull(triangle3D)


    # print "is in convex ?   {0}".format((not a.IsPntInConvexHull(non_hull,np.array([3.9,0,2.4]))) & a.IsPntInConvexHull(hull,np.array([0.2,2,2])))
    x_pnts = [1, 1.5, 2]
    y_pnts = [1, 1, 1]
    expected = [np.array([1.5, 1])]

    # a.PlotCvxHull(hull, cvx2=non_hull)
    a.ScatterPlot2D(x_pnts,
                    y_pnts,
                    "test",
                    ["x[m]", "y[m]"],
                    [(0,3.9),(0,4.07)],
                    cvx1=hull2d,
                    cvx2=non_hull2d,
                    expected=expected)
    err = []
    for i in range(len(x_pnts)):
        err.append(a.CalcErrorFromExpected([x_pnts[i], y_pnts[i]],expected))


    print "max error = {0}\nmin error = {1}\naverage error = {2}".format(max(err),min(err),np.average(err))

    # print hull.points
    # from src.wave2toa import Speaker
    # from src.wave2toa import recwav
    # rec = recwav()
    # rec.change_path('../inputs/blackmanharris5ms/1.WAV','in')
    # rec.PlotSignal('blackmanharris5ms.wav')
    # s = Speaker()
    # s.Define_ID(1)
    # s.BuildChirp()
    # sp_sig = {'signal': s.chirp, 'Fs': 88200, 'low_freq': 27500, 'high_freq': 30000, 'mat_signal': s.matlab_chirp}
    # rec_sig = {'signal': rec.signal, 'Fs': 88200, 'low_freq': 7500, 'high_freq': 30000, 'time': rec.rec_time}
    # b = SignalHandler()
    # c = SignalHandler()
    # b.defineParams(sp_sig)
    # c.defineParams(rec_sig)
    #
    # b.BPF(139)
    # c.BPF(139)
    # b.SmoothSig(139, step='a')
    # b.SmoothSig(139, step='b')
    # c.SmoothSig(139, step='a')
    # c.SmoothSig(139, step='b')
    # def CorrWith4(self, speakers):
    #     corr_list = []
    #     for speaker in speakers:
    #         corr_list.append(np.correlate(speaker.proccessed_signal.filtered_signal, speaker.chirp, 'full'))
    #     return corr_list
    # def FindTOA(self, corr_l, record, expected_signal_size, mode='max', thres=100):
    #     toa = []
    #     max_corr = []
    #     if mode == 'threshold':
    #         for v in corr_l:
    #             print "TBD " + str(thres)
    #     else:
    #         for v in corr_l:
    #             max_corr.append(max(v))
    #             max_index = v.argmax(axis=0)
    #             curr_toa = record.curr_time + (float(max_index) - expected_signal_size)/record.sample_rate
    #             toa.append(curr_toa)
    #
    #     return toa, max_corr

    print "Finish All"