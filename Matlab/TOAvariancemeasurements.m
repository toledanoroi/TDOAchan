% author: Roi Toledano
% this script checking the variance and expectation of time between TOA
% measurements. 
% Goal : 
%   reach variance < 50
% Improvement to do: 
%   1. look at 4-5 results and throw out the outliers
%   2. filtering the results with median/averaging filter size < 5 samples
%   3. record with higher pulse rate and merge each 5 samples to 1 sample.

clear all;
signal = audioread('/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/BH2.WAV');
expected_signals = load('/Users/roitoledano/git/TDOAchan/TDOAPostProcessing/inputs/allchirps.mat');
Fs = 500000;
tt = 1:1/Fs:5;
speed_of_voice = 343; % [m/s]
tested_sig = expected_signals.allchirp(:,1);
tested_sig = resample(tested_sig,1,2);

correlation = xcorr(signal,tested_sig);


figure; plot(correlation);

[pk,lk] = findpeaks(correlation,'MinPeakDistance',220000, 'MinPeakHeight' , 0.7 * max(correlation));

differ = diff(lk);
differ_amp = diff(pk);

mu = sum(differ)/length(differ);
a= differ - mu;
average_diff = sum(abs(a)) / length(a);
sigma = sqrt(sum(a.^2)/length(a));

delta_time = average_diff / Fs;
delta_dist = delta_time * speed_of_voice; % [s * m/s]

