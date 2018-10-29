import numpy as np
from scipy import signal
from matplotlib import pyplot as plt
from math import sin , cos , pi


def WaveformBuilder(params):
        tt = np.linspace(0, params['t1'], params['Fs'] * params['t1'])
        if params['type'] == 'builtin':
            tmp_chirp = signal.chirp(tt, params['f0'], params['t1'], params['f1'],
                                     method=params['mode'])
        elif params['type'] == 'unique':
            k = (params['f1'] - params['f0']) / (params['t1'])
            a = np.cos(tt)
            f = params['f0'] + k * tt + a
            tmp_chirp = np.sin(2 * pi * f * tt)
        if params['nrm']:
            nrm = np.linspace(1, 0.02, num=len(tmp_chirp))
            final_sig = np.multiply(nrm, tmp_chirp)
        else:
            final_sig = tmp_chirp

        return final_sig


if __name__ == '__main__':
    sig1_p = {'type': 'unique', 'mode':'linear', 'Fs': 192000, 't1': 0.1, 'nrm': False, 'f0': 30000, 'f1': 20000}
    sig2_p = {'type': 'builtin', 'mode': 'linear', 'Fs': 192000, 't1': 0.1, 'nrm': False, 'f0': 30000, 'f1': 20000}
    sig3_p = {'type': 'builtin', 'mode': 'hyperbolic', 'Fs': 192000, 't1': 0.1, 'nrm': False, 'f0': 30000, 'f1': 20000}
    sig4_p = {'type': 'builtin', 'mode': 'logarithmic', 'Fs': 192000, 't1': 0.1, 'nrm': False, 'f0': 30000, 'f1': 20000}
    sig5_p = {'type': 'builtin', 'mode': 'quadratic', 'Fs': 192000, 't1': 0.1, 'nrm': False, 'f0': 30000,
              'f1': 20000}


    list_params = [sig1_p,sig2_p,sig3_p,sig4_p,sig5_p]
    list_sig = []
    list_corr = []
    list_abgn = []
    for param in list_params:
        tmp = WaveformBuilder(param)
        list_sig.append(tmp)
        list_corr.append(np.correlate(tmp,tmp,'full'))
    t = np.linspace(0,0.02,len(list_corr[0]))


    fig = plt.figure()
    ax1 = fig.add_subplot(231)
    ax1.semilogy(np.linspace(0,0.2,len(list_corr[0])),list_corr[0],'r-')
    ax1.set_title("NLFM - cosine")
    ax2 = fig.add_subplot(232)
    ax2.semilogy(np.linspace(0,0.2,len(list_corr[1])),list_corr[1])
    ax2.set_title("linear")
    ax3 = fig.add_subplot(233)
    ax3.semilogy(np.linspace(0,0.2,len(list_corr[2])),list_corr[2])
    ax3.set_title("hyperbolic")
    ax4 = fig.add_subplot(234)
    ax4.semilogy(np.linspace(0,0.2,len(list_corr[3])),list_corr[3])
    ax4.set_title("logarithmic")
    ax5 = fig.add_subplot(235)
    ax5.semilogy(np.linspace(0,0.2,len(list_corr[4])),list_corr[4])
    ax5.set_title("Quadratic")

    plt.show()



    fig = plt.figure()
    ax1 = fig.add_subplot(231)
    ax1.plot(np.linspace(0,0.1,len(list_sig[0])),list_sig[0],'r-')
    ax1.set_title("NLFM - cosine")
    ax2 = fig.add_subplot(232)
    ax2.plot(np.linspace(0,0.1,len(list_sig[1])),list_sig[1])
    ax2.set_title("linear")
    ax3 = fig.add_subplot(233)
    ax3.plot(np.linspace(0,0.1,len(list_sig[2])),list_sig[2])
    ax3.set_title("hyperbolic")
    ax4 = fig.add_subplot(234)
    ax4.plot(np.linspace(0,0.1,len(list_sig[3])),list_sig[3])
    ax4.set_title("logarithmic")
    ax5 = fig.add_subplot(235)
    ax5.plot(np.linspace(0,0.1,len(list_sig[4])),list_sig[4])
    ax5.set_title("Quadratic")

    plt.show()