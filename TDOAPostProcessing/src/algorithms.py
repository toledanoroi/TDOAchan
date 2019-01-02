# project: TDOA based ultrasonic sound localization system
# authors: Roi Toledano & Yarden Avraham
# lab    : Bats lab
# guide  : PhD Yossi Yovel

import numpy as np
import scipy
from scipy import io
from scipy.spatial import distance
from time import time
from os.path import isfile
from termcolor import colored
import itertools
from collections import OrderedDict as OD
from matplotlib import pyplot as plt
from src.utils import UTILS

plotting = False
class ChanAlgo():

    def __init__(self, avg_dim=5, use_avg=False, Temprature_meas=''):
        self.wakeuptime = time()
        self.location_list = []

        self.dim = 3
        self.res = 0.01  # Tracking Resolution
        self.avg_dim = avg_dim
        self.use_avg = use_avg
        self.utils = UTILS()

        if isfile(Temprature_meas) & Temprature_meas.endswith('.mat'):
            temprature_mat = io.loadmat(Temprature_meas)
            mu_temprature = temprature_mat['mu']
            self.speedofvoice = 331.3 + 0.606 * mu_temprature
        else:
            self.speedofvoice = 343.21

    def DDOA(self, TimeToSpeaker, SpeakerLocations, alreadyTDOA=False):
        ddoa = scipy.array([])

        # A = np.append(SpeakerLocations, TimeToSpeaker, axis=1)
        #
        # TOA = np.zeros(4)
        if alreadyTDOA:
            Rnm = TimeToSpeaker
        else:
            Rnm = [(TimeToSpeaker[i] - TimeToSpeaker[0]) * self.speedofvoice for i in range(4)]
        # for i in range(0, 4):
        #     TOA[i] = TimeToSpeaker[i] - TimeToSpeaker[0]
        #     # TOA[i] = TimeToSpeaker[i] - min(TimeToSpeaker)
        #     Rnm.append(TOA[i] * self.speedofvoice)
        # # print TOA
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

    def FilterRoom(self, location_list, room_dimension):
        '''
        check if the algorithm result is out of room dimensions range and filter it.
        :param location_list: algorithm output : [<timestamp>,[<x>,<y>,<z>], <error>].
                              error used only for LUT room matrix algorithm
        :param room_dimension: list of 3 tuples, define the room size [(minx,maxx), (miny,maxy), (minz,maxz)]
        :return: loc_list_filtered - filtered results
        '''
        loc_list_filtered = []
        for loc in location_list:
            if (((loc[1][0] <= room_dimension[0][1]) & (loc[1][0] >= room_dimension[0][0])) &
                    ((loc[1][1] <= room_dimension[1][1]) & (loc[1][1] >= room_dimension[1][0])) &
                    ((loc[1][2] <= room_dimension[2][1]) & (loc[1][2] >= room_dimension[2][0]))):
                loc_list_filtered.append(loc)
            else:
                print "The result exceeds the room's dimensions\n\t {}".format(loc)
        return loc_list_filtered

    def chan_main(self, sp2mic, sp_location, timestamps, avg_list):
        '''
        main function of chan's algorithm . holds the algorithm sequence.
        :param sp2mic: list of speaker TOAs to the microphone. [[<toa1>,<toa2>,<toa3>,<toa4>],...]
        :param sp_location: speakers location in meters
        :param timestamps: list of TOA events timestamp  (-5ms from averaged toa)
        :param avg_list: int list, if needs averaging results. define the groups to average
                        [<g1>,<g2>,<g3>..]
        :return: list of all locations : [[<t1>,[<x>,<y>,<z>], <error>],...,[<t_final>,[<x>,<y>,<z>], <error>]].
                 error used only for LUT room matrix algorithm.
        '''

        # validity check

        cols = len(sp2mic[0])
        if len(sp_location) != 4:
            print colored("Error!! - Chan's Algorithm support only 4 speakers")
        # add here averaging results after throwing them.
        if self.use_avg:
            TDOA_dict = self.utils.Sp2MicToTDOA(sp2mic)
            measuredTDOA_vectors, timestamps = self.utils.AveragingSamples(avg_list, timestamps, self.avg_dim, TDOA_dict)
        for i in range(cols):
            if self.use_avg:
                curr_tdoa = [[measuredTDOA_vectors['TDOA_from_sp' + str(sp_number)][i][k]
                              for k in range(len(measuredTDOA_vectors['TDOA_from_sp1'][i]))]
                             for sp_number in range(1, 5)]
                D = self.DDOA(curr_tdoa[0], sp_location, alreadyTDOA=True)

            else:
                current_toa = [sp2mic[k][i] for k in range(len(sp2mic))]
                D = self.DDOA(current_toa,sp_location)

            POS = self.doChanForFour(D)
            print POS

            self.location_list.append([timestamps[i], POS[0], 0])

        return self.location_list

class RoomMatrix(object):
    def __init__(self, x, y, z, sp_list, avg_dim=5, res=0.1, Temprature_meas='', constant_z=-1):
        '''
        Initiate Room Matrix algorithm, define all parameters and create room quantization
        :param x: room length - x axis
        :param y: room length - y axis
        :param z: room length - z axis
        :param sp_list: list of all speakers objects
        :param avg_dim: groups size for averaging results
        :param res: resolution [m]  default = 0.1 [m] = 10 [cm]
        '''
        self.xlim = x
        self.ylim = y

        self.resolution = res
        self.speakers = {}
        self.EuclideanDistance = {}
        self.wakeup_time = time()
        self.finish_time = 0
        self.avg_dim = avg_dim
        self.utils = UTILS()

        self.dimx = int(np.ceil(self.xlim/self.resolution)) + 1
        self.dimy = int(np.ceil(self.ylim/self.resolution)) + 1
        if constant_z == -1:
            self.zlim = z
            self.dimz = int(np.ceil(self.zlim/self.resolution)) + 1
            z_t = np.linspace(0, self.zlim, self.dimz)
        else:
            self.zlim = constant_z
            self.dimz = 1
            z_t = np.array([constant_z])
        # define the room matrix
        x_t = np.linspace(0, self.xlim, self.dimx)
        y_t = np.linspace(0, self.ylim, self.dimy)

        # define room Quantization
        perm = np.array(list(itertools.product(x_t, y_t, z_t)))
        self.room_mat = perm.reshape(self.dimx, self.dimy, self.dimz, 3)

        #define speakers location
        for sp in sp_list:
            self.speakers[sp.id] = [sp.x, sp.y, sp.z]

        if isfile(Temprature_meas) & Temprature_meas.endswith('.mat'):
            temprature_mat = io.loadmat(Temprature_meas)
            mu_temprature = temprature_mat['mu']
            self.speedofvoice = 331.3 + 0.606 * mu_temprature
        else:
            self.speedofvoice = 343.21

    def CalcDistMatrix(self):
        '''
        This function calculates for each speaker the distance matrix
        from all the points in the room after quantization
        :return: generate self.mat4corr , expected distances matrix for correlation
        '''
        for key, value in self.speakers.items():
            tmp = np.tile(value, self.dimx * self.dimy * self.dimz)
            tmp2 = tmp.reshape(self.dimx, self.dimy, self.dimz, 3)
            # differ = (abs(tmp2 - self.room_mat)) ** 2
            #np.sqrt(differ[:, :, :, 0] + differ[:, :, :, 1] + differ[:, :, :, 2])
            self.EuclideanDistance['sp' + str(key)] = self.CalcEucDist2mat(tmp2, self.room_mat)
        shapeit = (self.EuclideanDistance['sp' + str(key)].reshape(self.dimx * self.dimy * self.dimz, 1)
                   for key in self.speakers.keys())
        self.mat4corr = np.column_stack(shapeit)

    def CalcTDOMat(self):
        self.TDOAmats = OD()
        tmp = self.mat4corr / float(self.speedofvoice)
        for i in range(len(self.speakers)):
            tmpvector = np.tile(tmp[:, i], len(self.speakers))
            tmpvector = np.transpose(tmpvector.reshape(len(self.speakers), len(tmp)))
            self.TDOAmats['TDOA_from_sp' + str(i+1)] = tmp - tmpvector

    def CalcEucDist2mat(self, mat1, mat2):
        differ = (abs(mat1 - mat2)) ** 2
        return np.sqrt(differ.sum(-1))

    def FindBestMatch(self, tdoa, mode='distance'):
        if mode == 'distance':
            tmp = np.tile(tdoa, self.dimx * self.dimy * self.dimz)
            tmp2 = tmp.reshape(self.dimx * self.dimy * self.dimz, len(self.speakers))
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

    def plotbytapsonly(self, sig, title):

        plt.plot(sig)
        plt.title(title)
        plt.grid()

    def WeightBestMatch(self, tdoa_list, mode='distance', consideration='all'):
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
                self.plotbytapsonly(dist[i], 'Room Euclidien Distance from measurement')
            legend.append(str(i))
        if plotting:
            plt.legend(legend)
            plt.show(block=False)

        # find shared best point:
        if consideration == 'all':
            summ = 0
            for di in dist:
                summ += di**2
            costfunction = np.sqrt(summ)
        elif consideration == '1':
            costfunction = dist[0]

        if plotting:
            plt.figure()
            self.plotbytapsonly(costfunction, 'Room Euclidien Distance from measurement')
            plt.show(block=False)

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

    def RoomMatMain(self, sp2mic, time_vect, avg_list, room_shape='square',use_avg=False, constant_z=-1):
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
        # if constant_z == -1:
        # calculate the Euclidean matrix. (LUT)
        self.CalcDistMatrix()
        # Generate relevant TDOAs matrices
        self.CalcTDOMat()
        # convert measured TOA to TDOA
        self.measuredTDOA_vectors = self.utils.Sp2MicToTDOA(sp2mic)
        # averaging if needed
        if use_avg:
            self.measuredTDOA_vectors, time_vect = self.utils.AveragingSamples(avg_list, time_vect, self.avg_dim, self.measuredTDOA_vectors)

        #create matrix for match filter
        for i in range(len(time_vect)):
            curr_tdoas = [self.measuredTDOA_vectors['TDOA_from_sp' + str(key)][i] for key in self.speakers.keys()]

            # find best match in LUT
            # if considered only one tdoa calculation
            mic_location, location_error = self.WeightBestMatch(curr_tdoas)#, consideration='1')
            locations_list.append([time_vect[i], mic_location, location_error])

        self.finish_time = time()
        print "algorithm time : {}".format(self.finish_time - self.wakeup_time)

        return locations_list



