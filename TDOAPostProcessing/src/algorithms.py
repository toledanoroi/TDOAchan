import numpy as np
import scipy
from scipy.spatial import distance
from time import time
from termcolor import colored
import itertools


class ChanAlgo():

    def __init__(self):
        self.wakeuptime = time()
        self.location_list = []
        self.speedofvoice = 343.0
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
            print colored(a,'red')
            D = self.DDOA(a,sp_location)
            POS = self.doChanForFour(D)
            print colored(POS,'red')
            self.location_list.append(POS[0])

        return self.location_list

class TaylorLS():
    def __init__(self):
        self.wakeuptime = time()
        self.location_list = []
        self.speedofvoice = 343.0
        self.dim = 3
        self.res = 0.01  # Tracking Resolution

class RoomMatrix(object):
    def __init__(self):
        self.xlim = 0
        self.ylim = 0
        self.zlim = 0
        self.resolution = 0.01
        self.speakers = {}
        self.EuclideanDistance = {}
        self.speedofvoice = 343.0
        self.wakeup_time = time()
        self.finish_time = 0

    def DefineRoomSize(self,x,y,z,res):
        self.xlim = x
        self.ylim = y
        self.zlim = z
        self.resolution = res


        self.dimx = int(np.ceil(self.xlim/self.resolution))
        self.dimy = int(np.ceil(self.ylim/self.resolution))
        self.dimz = int(np.ceil(self.zlim/self.resolution))
        # define the room matrix
        x_t = np.linspace(0, self.xlim, self.dimx)
        y_t = np.linspace(0, self.ylim, self.dimy)
        z_t = np.linspace(0, self.zlim, self.dimz)

        # define room Quantization
        perm = np.array(list(itertools.product(x_t, y_t, z_t)))
        self.room_mat = perm.reshape(self.dimx, self.dimy, self.dimz, 3)



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

    def CalcEucDist2mat(self,mat1, mat2):
        differ = (abs(mat1 - mat2)) ** 2
        return np.sqrt(differ.sum(-1))

    def DefineSpeakers(self,sp_list):
        for sp in sp_list:
            self.speakers[sp.id] = [sp.x, sp.y, sp.z]

    def FindBestMatch(self, toa, mode='distance'):
        if mode == 'distance':
            tmp = np.tile(toa, self.dimx * self.dimy * self.dimz)
            tmp2 = tmp.reshape(self.dimx * self.dimy * self.dimz, 4)
            dist = self.CalcEucDist2mat(tmp2, self.mat4corr)
            best = np.argmin(dist)
        elif mode == 'corr':
            best = np.argmax(np.dot(self.mat4corr, toa))

        return self.indextolocation(best)


    def indextolocation(self,index):
        x = int(np.floor(index/(self.dimy * self.dimz)))
        y = int(np.mod(index, x)/self.dimz)
        z = int(np.mod(np.mod(index, x), y))
        return self.room_mat[x, y, z, :]


    def RoomMatMain(self,sp2mic, speakers,room_size, resolution, room_shape='square'):
        '''

        :param sp2mic: TOA samples from each speaker to the microphone
        :param speakers: list of 4 speakers. wave2toa.Speaker()
        :param room_size: dictionary of the sizes in each axis of the room. {'x': ,'y':, 'z':} squaring the room
        :param resolution: quantization quality
        :param room_shape:  if not square, canceling part of the defined room.
        :return: locations list of the microphone.
        '''
        cols = len(sp2mic[0])
        timer = 0
        locations_list = []
        # initiate the room
        self.DefineRoomSize(room_size['x'], room_size['y'], room_size['z'], resolution)
        # self.DefineShape(room_shape)  # need to write if the room is not square.
        self.DefineSpeakers(speakers)
        # calculate the Euclidean matrix. (LUT)
        self.CalcDistMatrix()
        #create matrix for match filter

        for timestamp in range(cols):
            # need to create sp2mic - relevant timestamp   [TBD]

            current_toa = np.array([sp2mic[0][timestamp], sp2mic[1][timestamp],
                                    sp2mic[2][timestamp], sp2mic[3][timestamp]])*self.speedofvoice
            print colored(current_toa, 'red')

            # find best match in LUT
            mic_location = self.FindBestMatch(current_toa)
            locations_list.append([timer, mic_location])

        self.finish_time = time()

        return locations_list



if __name__ == '__main__':

    import wave2toa as wt
    a = RoomMatrix()
    a.DefineRoomSize(1,1,1,0.25)
    speakers = [wt.Speaker(), wt.Speaker(), wt.Speaker(), wt.Speaker()]
    for i in range(len(speakers)):
        speakers[i].Define_ID(i+1)

    # speakers[0].DefineLocation('s',0.0,0.0,1.0)
    # speakers[1].DefineLocation('s', 0.56346, 3.1346, 1.35)
    # speakers[2].DefineLocation('s', 2.1560, 0.3170, 1.370)
    # speakers[3].DefineLocation('s', 2.10, 3.15730, 1.45)

    speakers[0].DefineLocation('s',0.0,0.0,1.0)
    speakers[1].DefineLocation('s', 0.5, 3.0, 1.0)
    speakers[2].DefineLocation('s', 2.0, 0.0, 1.0)
    speakers[3].DefineLocation('s', 2.0, 3.0, 1.5)


    a.DefineSpeakers(speakers)
    a.CalcDistMatrix()
    # for key, item in a.EuclideanDistance.items():
    #     print key
    #     print item

    b = a.EuclideanDistance['sp1'].reshape(a.dimx * a.dimy * a.dimz, 1)
    c = a.EuclideanDistance['sp2'].reshape(a.dimx * a.dimy * a.dimz, 1)
    d = a.EuclideanDistance['sp3'].reshape(a.dimx * a.dimy * a.dimz, 1)
    e = a.EuclideanDistance['sp4'].reshape(a.dimx * a.dimy * a.dimz, 1)
    roi = np.column_stack((a.EuclideanDistance['sp1'].reshape(a.dimx * a.dimy * a.dimz, 1),
                           a.EuclideanDistance['sp2'].reshape(a.dimx * a.dimy * a.dimz, 1),
                           a.EuclideanDistance['sp3'].reshape(a.dimx * a.dimy * a.dimz, 1),
                           a.EuclideanDistance['sp4'].reshape(a.dimx * a.dimy * a.dimz, 1)))
    print a.EuclideanDistance['sp1']
    print b, c, d, e, roi
