# project: TDOA based ultrasonic sound localization system
# authors: Roi Toledano & Yarden Avraham
# lab    : Bats lab
# guide  : PhD Yossi Yovel

from src.wave2toa import RxMain
import numpy as np
from termcolor import colored

def ToaGenerator(base_path, params, prefix='5mili'):
    '''
    runs the algorithm for for all WAV files in a folder with required prefix
    generates TOA files for each wav
    :param base_path: the folder we want to search in.
    :param params: All tests params to send to RxMain function (imported from src.wave2toa)
    :param prefix: The prefix for wanted files
    '''
    import fnmatch
    import os
    err2d = []
    err3d = []
    matches = []
    # a = os.walk('/Users/roitoledano/Desktop/WLAN_vids')
    # for dirs in a:
    #     print dirs
    for root, dirnames, filenames in os.walk(base_path):
        for filename in fnmatch.filter(filenames, prefix + '*.WAV'):
            matches.append(os.path.join(root, filename))
    params['only_toa'] = True
    for match in matches:
        params['record_path'] = match
        RxMain(params)

def RunMultipleRecords(fin,params):
    '''
    take record list from the csv described in fin and runs test for each line.
    :param fin: the records list csv file path
    :param params: parameters of the test ,
    this function changes them each iteration and send it to RxMain function
    :return: dictionaries of errors (2D and 3D) from expected for each record
    avgerr["2D"], avgerr["3D"]
    '''
    from pandas import read_csv
    import os.path as op

    resfile = read_csv(fin)

    run_list = [(resfile['record_path'][i], resfile['matlab_path'][i], resfile['e_x'][i], resfile['e_y'][i], resfile['e_z'][i])
            for i in range(len(resfile['toa_path'])) if op.isfile(resfile['toa_path'][i])]
    avgerr = {"2D": {}, "3D": {}}
    for i in range(len(run_list)):
        params['matlab_path'] = run_list[i][1]
        params['record_path'] = run_list[i][0]
        params['expected_points'] = [np.array([run_list[i][2], run_list[i][3], run_list[i][4]])]
        params['expected_points2d'] = []
        for pnt in params['expected_points']:
            params['expected_points2d'].append(pnt[0:2])
        avgerr["2D"][params['record_path']], avgerr["3D"][params['record_path']] = RxMain(params)
    return avgerr["2D"], avgerr["3D"]

def ReadToaandExpected(fin):
    '''
    Takes TOA list from the csv described in fin and and returns its parameters.
    :param fin: the TOAs list csv file path
    :return: list of tuples [(toa_path, E[x], E[y], E[z]),...]
    '''
    from pandas import read_csv
    import os.path as op

    resfile = read_csv(fin)
    return [(resfile['toa_path'][i],resfile['e_x'][i],resfile['e_y'][i],resfile['e_z'][i])
            for i in range(len(resfile['toa_path'])) if op.isfile(resfile['toa_path'][i])]

def RunToaList(fin, params):
    '''
    Takes TOA list from the csv described in fin and run a test for each parameters line.
    :param fin: the TOAs list csv file path
    :param params: parameters of the test ,
    this function changes them each iteration and send it to RxMain function
    :return: dictionaries of errors (2D and 3D) from expected for each record
    err2d, err3d
    '''
    run_list = ReadToaandExpected(fin)
    for path, e_x, e_y, e_z in run_list:
        params['expected_points'] = [np.array([e_x, e_y, e_z])]
        params['expected_points2d'] = []
        for pnt in params['expected_points']:
            params['expected_points2d'].append(pnt[0:2])
        params['TOA_path'] = path
        print colored(path, 'red')
        err2d[path] = []
        err3d[path] = []
        for k in delt:
            params['resolution'] = k
            print colored("LUT room matrix Algorithm\n\tresolution = {0}\n\texpected"
                          " location = {1}".format(params['resolution'], params['expected_points'][0]), "red")

            a, b = RxMain(params)
            err2d[path].append(a)
            err3d[path].append(b)
    return err2d, err3d

def PlotErrorCurves(err2d, err3d, resolutions):
    '''
    plots The Euclidean errors vs algorithm resolution in meters
    :param err2d: dictionary of Euclidean errors 2D
    :param err3d: dictionary of Euclidean errors 3D
    :param resolutions: resolution list default = 0.01:0.2 jumps of 0.01 [m]
    '''
    from matplotlib import pyplot as plt
    plt.figure()
    txt = []
    for key, errors in err2d.items():
        plt.plot(resolutions, errors)
        plt.scatter(resolutions, errors)
        txt.append(key[key.rfind("/") + 1:])
    plt.legend(txt)
    plt.xlabel("resolution [m]")
    plt.ylabel("Euclidean error 2D [m]")
    plt.grid()
    plt.xticks(resolutions)
    plt.title("LUT room Matrix Algorithm Error vs matrix resolution - 2D")
    plt.show(block=False)
    plt.figure()
    txt = []
    for key, errors in err3d.items():
        plt.plot(resolutions, errors)
        plt.scatter(resolutions, errors)
        txt.append(key[key.rfind("/") + 1:])
    plt.legend(txt)
    plt.xlabel("resolution [m]")
    plt.ylabel("Euclidean error 3D [m]")
    plt.grid()
    plt.xticks(resolutions)
    plt.title("LUT room Matrix Algorithm Error vs matrix resolution - 3D")
    plt.show(block=False)
    mean_err2d = []
    mean_err3d = []
    plt.figure()
    kk = err2d.keys()
    for i in range(len(err2d[kk[0]])):
        mean_err2d.append(np.average([error2d[i] for error2d in err2d.values()]))
        mean_err3d.append(np.average([error3d[i] for error3d in err3d.values()]))
    plt.scatter(resolutions, mean_err2d)
    plt.scatter(resolutions, mean_err3d)
    plt.plot(resolutions, mean_err2d)
    plt.plot(resolutions, mean_err3d)
    plt.xlabel("resolution [m]")
    plt.ylabel("Euclidean error [m]")
    plt.grid()
    plt.xticks(resolutions)
    plt.title("LUT room Matrix Algorithm Error vs matrix resolution - averaged results")
    plt.legend(["2D", "3D"])
    plt.show()

if __name__ == '__main__':
    '''
    parameters definition for rx signal processing ultrasonic sound localization system
    ---------------------------------------------------------------------------------------------
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
    :param frecords - file path of all records to run
    :type string , must be a csv file with same format as TDOAchan/TDOAPostProcessing/inputs/records.csv
    ---------------------------------------------------------------------------------------------
    :param ftoas - file path of all TOA csv files to run
    :type string , must be a csv file with same format as TDOAchan/TDOAPostProcessing/toa_to_save/toacsvs.csv
    ---------------------------------------------------------------------------------------------
    :param frecbase - the folder with all .wav files to generate TOA for them.
    :type string , must be a path for real folder.
    ---------------------------------------------------------------------------------------------
    
        mode     |    operation
    _____________|_______________________________________________________________________________
        1        |  only generate TOAs from folder 
        2        |  run multiple records (WAV file) full - wave parsing and algorithm
        3        |  use toa csvs to run only the algorithm
        4        |  one shot - take one .wav file and run it full -> wave parsing and algorithm
    ---------------------------------------------------------------------------------------------
    '''

    # define parameters
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    res_iteration = False
    mode = 4
    params = {}
    algorithm_d = {'chan': 1,
                   'taylor': 2,
                   'room': 3,
                   'both': 4}
    # if TOA already exist
    params['TOA_path'] = ''  # '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/output/toa_of_thegoodsig.csv' #''/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/output/toa_record_1543256067.csv'  # '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/output/toa_record_1542718040.csv'
    params['speakers_frequencies'] = {'1': [34000, 41000],
                                      '2': [42000, 49000],
                                      '3': [20000, 27000],
                                      '4': [27000, 34000]}
    '''
    {'1': [20000, 27000],
      '2': [27000, 34000],
      '3': [34000, 41000],
      '4': [42000, 49000]}
    '''
    params['speakers_locations_d'] = {'1': [0.07, 0.07, 2.30],
                                      '2': [3.92, 0.07, 2.23],
                                      '3': [0.75, 3.82, 2.15],
                                      '4': [3.92, 4.01, 1.49]}
    params['room_sizes'] = {'x': 3.99,
                            'y': 4.08,
                            'z': 2.48}
    params['only_toa'] = False

    params['chirp_time'] = 0.005
    params['filter_size'] = 1001
    params['matlab_path'] = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/mat_files/allchirps_ChangedFreq.mat'
    params['record_path'] = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/goodsig_254_267_089.WAV'
    params['signal_mode'] = 1
    frecords = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/toa_to_save/records.csv'
    ftoas = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/toa_to_save/toacsvs.csv'
    frecbase = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/New_Inputs'

    params['algorithm'] = algorithm_d['room']
    params['resolution'] = 0.05  # [m]

    params['use_averaging_before_calculation'] = True
    params['time_factor'] = 2
    params['avg_group_size'] = 5

    # room shape convex hull parameters
    params['room3D'] = np.array([np.array([0, 0, 0]),
                                 np.array([0, 0, 2.48]),
                                 np.array([3.99, 0, 0]),
                                 np.array([3.99, 0, 2.48]),
                                 np.array([3.99, 4.04, 0]),
                                 np.array([3.99, 4.08, 2.48]),
                                 np.array([0.68, 4.08, 0]),
                                 np.array([0.68, 4.08, 2.48]),
                                 np.array([0.68, 1.82, 0]),
                                 np.array([0.68, 1.82, 2.48]),
                                 np.array([0, 1.82, 0]),
                                 np.array([0, 1.82, 2.48]),
                                 ])
    params['triangle3D'] = np.array([np.array([0.68, 4.08, 0]),
                                     np.array([0.68, 4.08, 2.48]),
                                     np.array([0.68, 1.82, 0]),
                                     np.array([0.68, 1.82, 2.48]),
                                     np.array([0, 1.82, 0]),
                                     np.array([0, 1.82, 2.48]),
                                     ])
    params['room2D'] = []
    params['triangle2D'] = []
    for pnt in params['room3D']:
        params['room2D'].append(pnt[0:2])
    for pnt in params['triangle3D']:
        params['triangle2D'].append(pnt[0:2])

    # expected_point for the relevant test
    params['expected_points'] = [np.array([2.54, 2.67, 0.89])]
    params['expected_points2d'] = []
    for pnt in params['expected_points']:
        params['expected_points2d'].append(pnt[0:2])

    err2d = {}
    err3d = {}
    delt = np.linspace(0.01, 0.2, 20)

    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # run the relevant mode
    if mode == 4:
        if res_iteration:
            for k in delt:
                params['resolution'] = k
                print colored("LUT room matrix Algorithm\n\tresolution = {0}\n\texpected"
                              " location = {1}".format(params['resolution'], params['expected_points'][0]), "red")
                a, b = RxMain(params)
                err2d[params['record_path']].append(a)
                err3d[params['record_path']].append(b)
        else:
            a, b = RxMain(params)
    elif mode == 1:
        # only generate TOAs
        ToaGenerator(frecbase, params, prefix='5mili')
    elif mode == 2:
        #  run multiple records (WAV file) full - wave parsing and algorithm
        _2Derr, _3Derr = RunMultipleRecords(frecords, params)
    elif mode == 3:
        # run loop with many records to examine the algorithm performence
        err2d, err3d = RunToaList(ftoas, params)

    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # plotting if needed
    if res_iteration | (mode == 3):
        PlotErrorCurves(err2d, err3d, delt)
