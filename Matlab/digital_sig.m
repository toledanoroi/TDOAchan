function [ dig_sig ] = digital_sig( t, fc , f_delta, samples_per_digital, mode )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% check for the first mod(t,samples) ==0
rem = mod(length(t),samples_per_digital);
l1 = length(t) - rem;
N = round(l1 / samples_per_digital);
% generate binary random signal
if mode==1
    bin_seq = randi([0 1], N , 1);
else
    goldseq = comm.GoldSequence('FirstPolynomial','x^5+x^2+1',...
    'SecondPolynomial','x^5+x^4+x^3+x^2+1',...
    'FirstInitialConditions',[0 0 0 0 1],...
    'SecondInitialConditions',[0 0 0 0 1],...
    'Index',0,'SamplesPerFrame',N);
    bin_seq = goldseq();
end

for i=1:N
    if bin_seq(i)==0
        tmp_sig(1+(i-1)*samples_per_digital:(i)*samples_per_digital) = cos(2*pi*t(1+(i-1)*samples_per_digital:(i)*samples_per_digital)*(fc - f_delta));
    else
        tmp_sig(1+(i-1)*samples_per_digital:(i)*samples_per_digital) = cos(2*pi*t(1+(i-1)*samples_per_digital:(i)*samples_per_digital)*(fc + f_delta));
    end
end
% "zero  modulated padding"
add_rem = cos(2*pi*t(l1:length(t) - 1)*(fc - f_delta));
dig_sig = [tmp_sig, add_rem];
        
end

