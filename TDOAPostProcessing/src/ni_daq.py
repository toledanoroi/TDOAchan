import nidaqmx as daq
from time import sleep


def ConnectDaq():
    pass

def TxSignal(signal , device_name, ports, round_time):
    with daq.Task() as task:
        task.ao_channels.add_ao_voltage_chan(device_name + "/" + ports)
    task.write(signal, auto_start=True)
    sleep(0.1)


def RxSignal(signal , device_name, ports, round_time):
    pass