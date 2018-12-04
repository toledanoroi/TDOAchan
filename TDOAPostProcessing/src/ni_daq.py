# project: TDOA based ultrasonic sound localization system
# authors: Roi Toledano & Yarden Avraham
# lab    : Bats lab
# guide  : PhD Yossi Yovel

# we use matlab for it but we create all the functions needed in python.

import nidaqmx as daq
from time import sleep


def TxSignal(signal , device_name, ports, round_time):
    with daq.Task() as task:
        task.ao_channels.add_ao_voltage_chan(device_name + "/" + ports)
    task.write(signal, auto_start=True)
    sleep(round_time)


def RxSignal(signal, device_name, ports, round_time):
    pass

def TempratureMeasurement(device_name, ports, round_time, rounds):
    import matplotlib.pyplot as plt
    plt.ion()
    i = 0
    with daq.Task() as task:
        task.ai_channels.add_ai_thrmcpl_chan(device_name + '/ai' + ports)
        while i < rounds:
            data = task.read(number_of_samples_per_channel=1)
            plt.scatter(i.data[0], c='r')
            plt.scatter(i.data[1], c='b')
            plt.pause(round_time)
            i += 1
            print data


