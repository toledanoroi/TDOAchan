import numpy as np
import scipy
from scipy.spatial import distance
from time import time
from termcolor import colored
import itertools
from collections import OrderedDict as OD
from matplotlib import pyplot as plt

plotting = False
class ChanAlgo():

    def __init__(self):
        self.wakeuptime = time()
        self.location_list = []
        self.speedofvoice = 343.21
        self.dim = 3
        self.res = 0.01  # Tracking Resolution

    def DDOA(self, TimeToSpeaker, SpeakerLocations):
        tmp = []
        toa = {}
        ddoa = scipy.array([])

        # A = np.append(SpeakerLocations, TimeToSpeaker, axis=1)

        TOA = np.zeros(4)
        Rnm = []
        for i in range(0, 4):
            TOA[i] = TimeToSpeaker[i] - TimeToSpeaker[0]
            # TOA[i] = TimeToSpeaker[i] - min(TimeToSpeaker)
            Rnm.append(TOA[i] * self.speedofvoice)
        # print TOA
        Rnm = np.matrix(Rnm)
        ddoa = np.append(SpeakerLocations, Rnm.T, axis=1)
        return ddoa

    def doChanForFour(self, D):
        M = D[:, :3]
        P = M[1:] - M[0]
        A = -P.I
        # A = np.linalg.pinv(-P)

        R = D[1:, 3]
        R_squared = scipy.matrix(R.A * R.A)

        K = scipy.zeros((4, 1))
        for n in range(4):
            K[n] = M[n] * M[n].T

        B = (R_squared - K[1:] + K[0]) / 2
        E = A * R
        F = A * B

        a = 1 - (E.T * E).item()
        b = 2 * (M[0] * E - F.T * E).item()
        c = (2 * (M[0] * F) - F.T * F - K[0]).item()
        discr = b * b - 4 * a * c

        if (discr < 0):
            return []
        h = scipy.sqrt(discr).item()

        poslist = []

        A1 = h
        A2 = -h

        for i in (A1, A2):
            R0 = (i - b) / (2 * a)

            if (R0 >= 0):
                T = E * R0 + F
                # poslist.append((T.A.squeeze(), R0))

                poslist.append((T.A.squeeze()))
        return poslist

    def chan_main(self,sp2mic,sp_location):

        rows = len(sp2mic)
        cols = len(sp2mic[0])
        for i in range(cols):
            a = []
            a.append(sp2mic[0][i])
            a.append(sp2mic[1][i])
            a.append(sp2mic[2][i])
            a.append(sp2mic[3][i])
            print (a,'red')
            D = self.DDOA(a,sp_location)
            POS = self.doChanForFour(D)
            print (POS,'red')
            self.location_list.append(POS[0])

        return self.location_list

class TaylorLS():
    def __init__(self):
        self.wakeuptime = time()
        self.location_list = []
        self.speedofvoice = 343.21
        self.dim = 3
        self.res = 0.01  # Tracking Resolution

class RoomMatrix(object):
    def __init__(self,x,y,z,sp_list,avg_dim=5,res=0.1):
        self.xlim = x
        self.ylim = y
        self.zlim = z
        self.resolution = res
        self.speakers = {}
        self.EuclideanDistance = {}
        self.speedofvoice = 343.21
        self.wakeup_time = time()
        self.finish_time = 0
        self.avg_dim = avg_dim

        self.dimx = int(np.ceil(self.xlim/self.resolution)) + 1
        self.dimy = int(np.ceil(self.ylim/self.resolution)) + 1
        self.dimz = int(np.ceil(self.zlim/self.resolution)) + 1
        # define the room matrix
        x_t = np.linspace(0, self.xlim, self.dimx)
        y_t = np.linspace(0, self.ylim, self.dimy)
        z_t = np.linspace(0, self.zlim, self.dimz)
        # define room Quantization
        perm = np.array(list(itertools.product(x_t, y_t, z_t)))
        self.room_mat = perm.reshape(self.dimx, self.dimy, self.dimz, 3)

        #define speakers location
        for sp in sp_list:
            self.speakers[sp.id] = [sp.x, sp.y, sp.z]

    def CalcDistMatrix(self):
        for key, value in self.speakers.items():
            tmp = np.tile(value, self.dimx * self.dimy * self.dimz)
            tmp2 = tmp.reshape(self.dimx, self.dimy, self.dimz, 3)
            # differ = (abs(tmp2 - self.room_mat)) ** 2
            #np.sqrt(differ[:, :, :, 0] + differ[:, :, :, 1] + differ[:, :, :, 2])
            self.EuclideanDistance['sp' + str(key)] = self.CalcEucDist2mat(tmp2, self.room_mat)

        self.mat4corr = np.column_stack((self.EuclideanDistance['sp1'].reshape(self.dimx * self.dimy * self.dimz, 1),
                                         self.EuclideanDistance['sp2'].reshape(self.dimx * self.dimy * self.dimz, 1),
                                         self.EuclideanDistance['sp3'].reshape(self.dimx * self.dimy * self.dimz, 1),
                                         self.EuclideanDistance['sp4'].reshape(self.dimx * self.dimy * self.dimz, 1)))

    def CalcTDOMat(self):
        self.TDOAmats = OD()
        tmp = self.mat4corr / float(self.speedofvoice)
        for i in range(4):
            tmpvector = np.tile(tmp[:,i], 4)
            tmpvector = np.transpose(tmpvector.reshape(4,len(tmp)))
            self.TDOAmats['TDOA_from_sp' + str(i+1)] = tmp - tmpvector

    def CalcEucDist2mat(self,mat1, mat2):
        differ = (abs(mat1 - mat2)) ** 2
        return np.sqrt(differ.sum(-1))

    def FindBestMatch(self, tdoa, mode='distance'):
        if mode == 'distance':
            tmp = np.tile(tdoa, self.dimx * self.dimy * self.dimz)
            tmp2 = tmp.reshape(self.dimx * self.dimy * self.dimz, 4)
            dist = self.CalcEucDist2mat(tmp2, self.TDOAmats['TDOA_from_sp1'])
            if plotting:
                self.plotbytapsonly(dist, 'Room Euclidien Distance from measurement')
            best = np.argmin(dist)
            error = min(dist)
        elif mode == 'corr':
            k = np.dot(self.TDOAmats['TDOA_from_sp1'], tdoa)
            if plotting:
                self.plotbytapsonly(k, 'Room Euclidien Distance correlation from measurement')
            best = np.argmax(k)
            error = max(k)

        return self.indextolocation(best)

    def plotbytapsonly(self,sig, title):

        plt.plot(sig)
        plt.title(title)
        plt.grid()

    def WeightBestMatch(self, tdoa_list, mode='distance',consideration='all'):
        dist = []
        legend = []
        if plotting:
            plt.figure()
        for i in range(len(tdoa_list)):
            tmp = np.tile(tdoa_list[i], self.dimx * self.dimy * self.dimz)
            tmp2 = tmp.reshape(self.dimx * self.dimy * self.dimz, 4)
            Tdoa_mat_keys = self.TDOAmats.keys()
            if mode == 'distance':
                dist.append(self.CalcEucDist2mat(tmp2, self.TDOAmats[Tdoa_mat_keys[i]]))
            elif mode == 'corr':
                dist.append(np.dot(self.TDOAmats[Tdoa_mat_keys[i]], tdoa_list[i]))
            if plotting:
                self.plotbytapsonly(dist[i],'Room Euclidien Distance from measurement')
            legend.append(str(i))
        if plotting:
            plt.legend(legend)
            plt.show()

        # find shared best point:
        if consideration == 'all':
            costfunction = np.sqrt(dist[0]**2 + dist[1]**2 + dist[2]**2 + dist[3]**2)
        elif consideration == '1':
            costfunction = dist[0]

        if plotting:
            plt.figure()
            self.plotbytapsonly(costfunction, 'Room Euclidien Distance from measurement')
            plt.show()

        if mode == 'distance':
            best = np.argmin(costfunction)
        elif mode == 'corr':
            best = np.argmax(costfunction)
        error = costfunction[best]

        return self.indextolocation(best), error

    def indextolocation(self,index):
        x = int(np.floor(index/(self.dimy * self.dimz)))
        y = int(np.floor((index - x * (self.dimy * self.dimz))/self.dimz))
        z = int(index - x * (self.dimy * self.dimz) - y * self.dimz)
        return self.room_mat[x, y, z, :]

    def Sp2MicToTDOA(self,sp2mic,avg_list,use_avg,time_vect):
        rows = len(sp2mic)
        self.measuredTDOA_vectors = OD()
        for i in range(rows):
            tmpvector = np.transpose(np.array(sp2mic) - np.array(sp2mic[i]))
            self.measuredTDOA_vectors['TDOA_from_sp' + str(i+1)] = tmpvector
        if use_avg:
            self.measuredTDOA_vectors, time_vect = self.AveragingSamples(avg_list,time_vect)
        return time_vect

    def AveragingSamples(self,avg_list,time_vect):
        from statistics import mean
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
            curr_ += self.avg_dim
            for i in range(len(v)):
                v[i] = self.measuredTDOA_vectors['TDOA_from_sp' + str(i+1)][last:curr]
                v_avg['TDOA_from_sp' + str(i+1)].append([mean([k[j] for k in v[i]]) for j in range(len(v[i][0]))])
            time_avg.append(mean(time_vect[last_:curr_]))
            last = curr
            last_ = curr_

        return v_avg, time_avg

    def RoomMatMain(self, sp2mic, time_vect, avg_list, room_shape='square',use_avg=False):
        '''
        :param sp2mic: TOA samples from each speaker to the microphone
        :param time_vect:
        :param avg_list:
        :param room_shape:  if not square, canceling part of the defined room.
        :param use_avg:
        :return: locations list of the microphone.
        '''

        cols = len(sp2mic[0])
        locations_list = []

        # calculate the Euclidean matrix. (LUT)
        self.CalcDistMatrix()
        # Generate relevant TDOAs matrices
        self.CalcTDOMat()
        # convert measured TOA to TDOA and averaging if needed
        time_vect = self.Sp2MicToTDOA(sp2mic,avg_list, use_avg,time_vect)

        #create matrix for match filter

        for i in range(len(time_vect)):
            current_tdoa1 = self.measuredTDOA_vectors['TDOA_from_sp1'][i]
            current_tdoa2 = self.measuredTDOA_vectors['TDOA_from_sp2'][i]
            current_tdoa3 = self.measuredTDOA_vectors['TDOA_from_sp3'][i]
            current_tdoa4 = self.measuredTDOA_vectors['TDOA_from_sp4'][i]

            # find best match in LUT
            # if considered only one tdoa calculation
            mic_location, location_error = self.WeightBestMatch([current_tdoa1,current_tdoa2,current_tdoa3,current_tdoa4])#, consideration='1')
            locations_list.append([time_vect[i], mic_location, location_error])


        self.finish_time = time()
        print "algorithm time : {}".format(self.finish_time - self.wakeup_time)

        return locations_list



