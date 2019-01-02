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

    for root, dirnames, filenames in os.walk(base_path):
        for filename in fnmatch.filter(filenames, prefix + '*.wav'):
            match = os.path.join(root, filename)
            params['only_toa'] = True
            params['record_path'] = match
            params['point_name'] = filename[:filename.rfind('.')]
            RxMain(params)

def RunWithoutExpected(base_path, params, prefix='a',
                       locations_path='/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/output/location_results',
                       ToA_path = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/output/toas_samples'):
    '''
        runs the algorithm for for all WAV files in a folder with required prefix
        generates TOA files for each wav and generate locations list
        :param base_path: the folder we want to search in.
        :param params: All tests params to send to RxMain function (imported from src.wave2toa)
        :param prefix: The prefix for wanted files
        '''
    import fnmatch
    import os
    from src.utils import UTILS

    for path in [locations_path, ToA_path]:
        try:
            os.rmdir(path)
        except OSError:
            print colored("removing the directory %s failed" % path, 'blue')
        try:
            os.mkdir(path)
        except OSError:
            if os.path.isdir(path):
                for root, dirnames, filenames in os.walk(path):
                    for filename in fnmatch.filter(filenames, prefix + '*.wav'):
                        match = os.path.join(root, filename)
                        os.remove(match)
                print ("All the files in the directory %s removed" % path)
            else:
                print ("Creation of the directory %s failed" % path)
        else:
            print ("Successfully created the directory %s " % path)

    params['expected_points'] = [np.array([-1, -1, -1])]
    ut = UTILS()


    for root, dirnames, filenames in os.walk(base_path):
        for filename in fnmatch.filter(filenames, prefix + '*.wav'):
            match = os.path.join(root, filename)
            params['record_path'] = match
            params['point_name'] = filename[:filename.rfind('.')]
            params['ToAs_file'] = ToA_path + '/' + params['point_name'] + '.csv'
            params['unique_path'] = locations_path + '/location_res_' + params['point_name'] + '.csv'
            RxMain(params)
    output = locations_path + '/merged_results.csv'
    results = ut.MergecsvAndGenerateForPlotting(locations_path, prefix='location_',output_path=output)

    if params['constant_z'] != -1:
        PlotResults(output, params)
    else:
        PlotResults(output, params, show='3D')














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

def PlotResults(fin, _params, show='2D'):
    from src.utils import UTILS as ut
    from pandas import read_csv , DataFrame
    utils_obj = ut()
    # define room convex
    hull, non_hull, hull2d, non_hull2d = utils_obj.DefineRoom(_params['room3D'],
                                                              _params['triangle3D'],
                                                              shapewithnonedgepoint=False,
                                                              plotting=False
                                                             )
    # read file:
    results = read_csv(fin)
    if ((results['E[x]'][0] >= 0) & (results['E[y]'][0] >= 0) & (results['E[z]'][0] >= 0) ):
        expected_locations = [[results['E[x]'][i], results['E[y]'][i]]
                              for i in xrange(len(results['E[x]']))]
    else:
        expected_locations = None

    if show == '2D':
        utils_obj.ScatterPlot2D(results['X [m]'],
                                results['Y [m]'],
                                'Room LUT algorithm results - 2D',
                                ['X [m]', 'Y [m]'],
                                [(0 ,_params['room_sizes']['x']), (0, _params['room_sizes']['y'])],
                                cvx1=hull2d,
                                cvx2=non_hull2d,
                                expected=expected_locations,
                                points=results['point_set']
                                )
    else:
        utils_obj.ScatterPlot3D(results['X [m]'],
                                results['Y [m]'],
                                results['Z [m]'],
                                'Room LUT algorithm results - 3D',
                                ['X [m]', 'Y [m]', 'Z [m]'],
                                [(0, _params['room_sizes']['x']), (0, _params['room_sizes']['y']), (0, _params['room_sizes']['z'])],
                                cvx1=hull,
                                cvx2=non_hull,
                                expected=expected_locations,
                                points=results['point_set']
                                )

    # plotting error
    if expected_locations == None:
        from matplotlib.pyplot import show
        show()
    else:
        from matplotlib.pyplot import plot, show, xlabel, ylabel, title, scatter, grid
        plot(results['Iteration'], results['Error'])
        scatter(results['Iteration'], results['Error'])
        xlabel("Iteration number")
        ylabel("Euclidean Error")
        title("Euclidian distance error per iteration")
        grid()
        show()

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
    from matplotlib.pyplot import show
    show()
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
        5        |  plot results from results file.
        6        |  run .wav files from folder without expected result
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
    params['TOA_path'] = '/005.m'  # '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/output/toa_of_thegoodsig.csv' #''/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/output/toa_record_1543256067.csv'  # '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/output/toa_record_1542718040.csv'
    params['number_of_speakers'] = 4
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
    # params['speakers_locations_d'] = {'1': [0.323, 3.646, 0.504],
    #                                   '2': [1.75, 0.43, 0.589],
    #                                   '3': [1.738, 2.843, 0.909],
    #                                   '4': [0.059, 0.158, 2.21]}

    params['speakers_locations_d'] = {'1': [0.2973, 3.6627, 0.504],
                                      '2': [1.8253, 0.4246, 0.589],
                                      '3': [1.7077, 2.8073, 0.9979],
                                      '4': [0.1039, 0.2065, 2.21]}

    params['room_sizes'] = {'x': 1.966,
                            'y': 4.272,
                            'z': 3.051}

    params['constant_z'] = 0.61
    params['only_toa'] = False
    import time
    params['point_name'] = '../output/TOA_' + str(int(time.time())) + '.csv'

    params['chirp_time'] = 0.002
    params['filter_size'] = 1001
    params['matlab_path'] = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/chirp_2m_bh.mat'
    params['record_path'] = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/batroom/a0000002.wav'
    params['signal_mode'] = 1
    params['unique_path'] = 'no_path'
    params['ToAs_file'] = 'no_path'
    frecords = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/toa_to_save/records.csv'
    ftoas = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/toa_to_save/toacsvs_digital.csv'
    frecbase = '/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/chirp_2m_bh'

    params['algorithm'] = algorithm_d['room']
    params['resolution'] = 0.005  # [m]
    params['use_averaging_before_calculation'] = True
    params['time_factor'] = 2
    params['avg_group_size'] = 5

    # room shape convex hull parameters
    params['room3D'] = np.array([np.array([0, 0, 0]),
                                 np.array([0, 0, params['room_sizes']['z']]),
                                 np.array([params['room_sizes']['x'], 0, 0]),
                                 np.array([params['room_sizes']['x'], 0, params['room_sizes']['z']]),
                                 np.array([0, params['room_sizes']['y'], 0]),
                                 np.array([0, params['room_sizes']['y'], params['room_sizes']['z']]),
                                 np.array([params['room_sizes']['x'], params['room_sizes']['y'], 0]),
                                 np.array([params['room_sizes']['x'], params['room_sizes']['y'], params['room_sizes']['z']])
                                 ])
    params['triangle3D'] = np.array([np.array([0.676, params['room_sizes']['y'], 0]),
                                     np.array([0.676, params['room_sizes']['y'], params['room_sizes']['z']]),
                                     np.array([0.676, 1.82, 0]),
                                     np.array([0.676, 1.82, params['room_sizes']['z']]),
                                     np.array([0, 1.82, 0]),
                                     np.array([0, 1.82, params['room_sizes']['z']])
                                     ])
    params['room2D'] = []
    params['triangle2D'] = []
    for pnt in params['room3D']:
        params['room2D'].append(pnt[0:2])
    for pnt in params['triangle3D']:
        params['triangle2D'].append(pnt[0:2])

    # expected_point for the relevant test
    params['expected_points'] = [np.array([1, 2, 0.61])]
    params['expected_points2d'] = []
    for pnt in params['expected_points']:
        params['expected_points2d'].append(pnt[0:2])

    err2d = {}
    err3d = {}
    # delt = np.linspace(0.01, 0.2, 20)
    delt = [0.005]
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # run the relevant mode
    if mode == 4:
        params['point_name'] = params['record_path'][params['record_path'].rfind('/') + 1: params['record_path'].rfind('.')]
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
        ToaGenerator(frecbase, params, prefix='a')
    elif mode == 2:
        #  run multiple records (WAV file) full - wave parsing and algorithm
        _2Derr, _3Derr = RunMultipleRecords(frecords, params)
    elif mode == 3:
        # run loop with many records to examine the algorithm performence
        err2d, err3d = RunToaList(ftoas, params)
    elif mode == 5:
        PlotResults('/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/toa_to_save/locations_results_17122018.csv',
                    params)
    elif mode == 6:
        RunWithoutExpected(frecbase, params, prefix='a')



    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------
    # plotting if needed
    if res_iteration | (mode == 3):
        PlotErrorCurves(err2d, err3d, delt)
