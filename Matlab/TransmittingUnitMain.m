clc
format long;

%% Test Parameters 

Fs=500000;  % 250000
signal_time = 0.002; % seconds
window = 'blackman_harris';  % 'blackman_harris' , 'hamming' , 'linear_diagonal' , 'nothing'
signal_type = 'chirp';
welcome_msg = 'Hi everyone I am the Transmitting Unit in the localization system of Roy Toledano and Yarden Avraham';
loops = 1000;
pause_time = 0.2;
freqs_mat = [20000 27000; 27000 34000; 34000 41000; 41000 49000;];
digital_distance = 3000;
gain = 5;  % linear gain
mode = 3;  
%           Mode 0 - Welcome                                 %
%           Mode 1 - Debug                                   %
%           Mode 2 - One by one                              %
%           Mode 3 - All together                            %
%           Mode 4 - channel frequency response measdurement %             
%           Mode 5 - digital binary fm modulated signal      %
%% functions calling
t=0:1/Fs:signal_time;

[sig1, sig2, sig3, sig4] = SpeakersWaveformBuilder( mode, signal_type, window, freqs_mat, t, signal_time,Fs,digital_distance);
allchirp = [sig1,sig2,sig3,sig4];

DaqHandler( mode ,loops, Fs, allchirp, pause_time, gain )




