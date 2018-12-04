# project: TDOA based ultrasonic sound localization system
# authors: Roi Toledano & Yarden Avraham
# lab    : Bats lab
# guide  : PhD Yossi Yovel

# just for debug!!!!

import numpy
import pylab
import scipy.io.wavfile
import scipy.signal

import utils

def short_time_energy_and_zero_cross_rate():
  """Example: Computation of ST-ZCR and STE of a speech signal."""

  fs, x = scipy.io.wavfile.read('../inputs/second_test.wav')
  x= numpy.array(x, dtype=float)
  x = x.flatten()
  t = numpy.arange(len(x)) * (1.0 / fs)
# x=[]
 # for couple in x2:
  #  x.append(couple[0]);



  # Find the short time zero crossing rate.
  zc = utils.stzcr(x, scipy.signal.get_window("boxcar", 201))
  
  # Find the short time energy.
  e = utils.ste(x, scipy.signal.get_window("hamming", 201))
  
  # pylab.figure()
  # pylab.subplot(311)
  # pylab.plot(t, x)
  # pylab.title('Speech signal (TEST1.wav)')
  # pylab.subplot(312)
  # pylab.plot(t, zc, 'r', linewidth=2)
  # pylab.title('Short-time Zero Crossing Rate')
  # pylab.subplot(313)
  # pylab.plot(t, e, 'm', linewidth=2)
  # pylab.xlabel('t (s)')
  # pylab.title('Short-time Energy')
  #
  # pylab.figure()
  # pylab.plot(t, x / x.max(), label="TEST1.wav")
  # pylab.hold(True)
  # pylab.plot(t, zc / zc.max(), 'r', label="zero crossing rate")
  # pylab.plot(t, e / e.max(), 'm', label="energy")
  # pylab.legend()


def run_examples():
  short_time_energy_and_zero_cross_rate()
  pylab.show()

def run_delay():
  a = numpy.zeros(2000)
  a[50:1500] = 1

  h = scipy.signal.firwin(numtaps=10001, cutoff=[60000], nyq=125000)
  filtered_signal = scipy.signal.convolve(a, h)
  corr_sig = numpy.correlate(a,numpy.array([0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]),'full')
  corr_a = numpy.correlate(a, a, 'full')

  print len([0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0])

  scipy.signal.chirp()

  pylab.plot(a)
  # pylab.plot(filtered_signal)
  pylab.plot(corr_sig)
  pylab.show()

  for i in range(len(filtered_signal)):
    if (corr_sig[i] > 0.5):
      s = i

  print s




if __name__ == "__main__":
  # run_examples()
  run_delay()


