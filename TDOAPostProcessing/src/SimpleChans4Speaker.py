import scipy
import scipy.stats
import scipy.linalg
import numpy as np
import time
import csv
import wave

speedofsound = 343.0
res = 0.01  # Tracking Resolution
dim = 3  # Tracking Dimensions


def DDOA(TimeToSpeaker, SpeakerLocations):
    tmp = []
    toa = {}
    ddoa = scipy.array([])

    A = np.append(SpeakerLocations, TimeToSpeaker, axis=1)

    TOA = np.zeros(4)
    Rnm = []
    for i in range(0, 4):
        TOA[i] = TimeToSpeaker[i] - TimeToSpeaker[0]
        #TOA[i] = TimeToSpeaker[i] - min(TimeToSpeaker)
        Rnm.append(TOA[i] * speedofsound)
    #print TOA
    Rnm = np.matrix(Rnm)
    ddoa = np.append(SpeakerLocations, Rnm.T, axis=1)
    return (ddoa)


def doChanForFour(D):
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

    for i in (A1,A2):
        R0 = (i - b) / (2 * a)

        if (R0 >= 0):

            T = E * R0 + F
            #poslist.append((T.A.squeeze(), R0))

            poslist.append((T.A.squeeze()))
    return poslist


Speaker1_x = 0.0
Speaker1_y = 0.0
Speaker1_z = 0.0

Speaker2_x = 0.30475778  # -0.055 #random.uniform(6.8,7.1) # 7.047
Speaker2_y = 6.12861056
Speaker2_z = 0.24361689

Speaker3_x = 8.71575555  # random.uniform(6.8,7.1) #6.842 #7.087
Speaker3_y = 7.68148882  # random.uniform(3.8,4.2) #3.82 #4.119
Speaker3_z = 0.94786345

Speaker4_x = 8.89878992
Speaker4_y = 1.08542346  # random.uniform(3.8,4.2) #3.82 #4.2708
Speaker4_z = -0.66767797

with open('2017-06-09-15-36-27.csv') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    all_list = []

    for row in spamreader:
        temp_list = []
        row = row[0].split(',')
        # print row
        for i in range(0, 4):
            temp_list.append(float(row[i]))

        all_list.append(temp_list)

SpeakerLocations = np.matrix(
    [[Speaker1_x, Speaker1_y, Speaker1_z], [Speaker2_x, Speaker2_y, Speaker2_z], [Speaker3_x, Speaker3_y, Speaker3_z],
     [Speaker4_x, Speaker4_y, Speaker4_z]])
start = time.time()

Position_lists = []


for i in range(0, 57):
    start_time = time.time()

    TimeToSpeaker = np.array([[all_list[i][0]], [all_list[i][1]], [all_list[i][2]], [all_list[i][3]]])

    start = time.time()

    D = DDOA(TimeToSpeaker, SpeakerLocations)

    POS = doChanForFour(D)

    print POS

    # plt.plot(X,Y)
    # plt.show()
    try:
        end = time.time()

        Position_lists.append(POS[0])
        #print("--- %s seconds ---" % (time.time() - start_time))

        end = time.time()
        #print "Time " + str(end - start)
    except:
        pass
