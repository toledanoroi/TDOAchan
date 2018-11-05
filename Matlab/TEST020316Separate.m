%clear all;
clc

format long;

Fs=250000;
t=0:1/Fs:0.001;
%fileID = fopen('exptable.txt','w');
results = {};

locationresults = [];
timereults = {};

locationresultsofeachstep = [];
timeresultsforeachstep = {};

fid1 = fopen('../Logs/log.txt','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Define Mode                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Mode 0 - Welcome                   %
%           Mode 1 - Debug     
%
%           Mode 2 - One by one                %
%           Mode 3 - All together              %
Mode = 3;
gain = 5;
disp('The TX system working at mode:');
disp(Mode);

if Mode==1
    wav1=1*tts('One',1);
    wav2=1*tts('Two',1);
    wav3=1*tts('Three',1);
    wav4=1*tts('Four',1);
    wav1(numel(wav4)) = 0;
    wav2(numel(wav4)) = 0;
    wav3(numel(wav4)) = 0;
    max_s = max([length(wav1),length(wav2),length(wav3),length(wav4)]);
    wav1 =  [wav1; zeros(1,max_s - length(wav1)).'];
    wav2 =  [wav2; zeros(1,max_s - length(wav2)).'];
    wav3 =  [wav3; zeros(1,max_s - length(wav3)).'];
    wav4 =  [wav4; zeros(1,max_s - length(wav4)).'];
end

if Mode==0
    wav1= 5 * tts('Hi everyone I am the Transmitting Unit in the localization system of Roy Toledano and Yarden Avraham',1);
    wav2= wav1;
    wav3= wav1;
    wav4= wav1;
end

if Mode==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Speaker 1          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [data , fs_arian] = audioread('D:\RoiT\final_project\matlab\1.wav');    
    y1=data;
    y1 = chirp(t,300,1,5000);
    y1 = chirp(t,20000,1,30000);
    y1=y1'; %a chirp signal
%     nrm=linspace(1,0.02,length(y1)); %"convolving" the chirp for higher amplitudes in high frequencies
%     y1=3*y1.*nrm';
    disp('done')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Speaker 2          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y2=chirp(t,30000,0.001,35000);
    y2 = chirp(t,42000,0.001,37000);
    
    y2 = data;
    y2 = chirp(t,5000,1,10000);
    y2 = chirp(t,30000,1,40000);
    y2=y2'; %a chirp signal
%     nrm=linspace(1,0.02,length(y2)); %"convolving" the chirp for higher amplitudes in high frequencies
%     y2=3*y2.*nrm';
    disp('done')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Speaker 3          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y3=chirp(t,23000,0.001,28000);
    y3 = chirp(t,42000,0.001,37000);
    y3 = data;
    y3 = chirp(t,10000,1,15000);
    y3 = chirp(t,40000,1,50000);
    y3=y3'; %a chirp signal
%     nrm=linspace(1,0.02,length(y3)); %"convolving" the chirp for higher amplitudes in high frequencies
%     y3=3*y3.*nrm';
    disp('done')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Speaker 4          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y4= chirp(t,15000,1,20000);
    y4= chirp(t,50000,1,60000);
%     y4 = chirp(t,42000,0.001,37000);
%     y4 = data;
    y4=y4'; %a chirp signal
%     nrm=linspace(1,0.02,length(y4)); %"convolving" the chirp for higher amplitudes in high frequencies
%     y4=3*y4.*nrm';
    disp('done')
    
    allchirp = [y1 y2 y3 y4];
%     plotting for arian
%     figure; hold on;
%     plot(tttt,allchirp(1,:));
%     plot(tttt,allchirp(2,:));
%     plot(tttt,allchirp(3,:));
%     plot(tttt,allchirp(4,:));
    
    a = allchirp(:,1);
    b = allchirp(:,2);
    c = allchirp(:,3);
    d = allchirp(:,4);
    figure; hold on;
    plot(t,a);
    plot(t,b);
    plot(t,c);
    plot(t,d);
end

if Mode==3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Speaker 1          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y1=chirp(t,20000,0.001,28000);
    y11=chirp(t,28000,0.001,20000);
    y111=[y1,y11];
    y111=y111'; %a chirp signal
    %nrm=linspace(1,0.04,length(y111)); %"convolving" the chirp for higher amplitudes in high frequencies
    %y1fft=fft(y1,251);
    win1=hamming(length(y111));
    y111=y111.*win1;
	%y111=y111.*nrm';
    %y12=win1.*y1fft
    disp('done')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Speaker 2          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y2=chirp(t,32000,0.001,40000);
    y22=chirp(t,40000,0.001,32000);
    y222=[y2,y22];
    y222=y222'; %a chirp signal
    win2=hamming(length(y222));
    y222=y222.*win2;
    %nrm=linspace(1,0.04,length(y222)); %"convolving" the chirp for higher amplitudes in high frequencies
    %y222=y222.*nrm';
    disp('done')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Speaker 3          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y3= chirp(t,44000,0.001,52000);
    y33= chirp(t,52000,0.001,44000);
    y333=[y3,y33];
    y333=y333.';
    win3=hamming(length(y333));
    y333=y333.*win3;
    %nrm=linspace(1,0.04,length(y333)); %"convolving" the chirp for higher amplitudes in high frequencies
    %y333=y333.*nrm';
    disp('done')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%         Speaker 4          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y4=chirp(t,56000,0.001,64000);
    y44=chirp(t,64000,0.001,56000);
    y444= [y4,y44]
    y444=y444'; %a chirp signal
    win4=hamming(length(y444));
    y444=y444.*win4;
    %nrm=linspace(1,0.04,length(y444)); %"convolving" the chirp for higher amplitudes in high frequencies
	%y444=y444.*nrm';
    disp('done')
   
    % add zeros 
    zz = zeros(2500,1);
    y111 = [zz; y111; zz;];
    y222 = [zz; y222; zz;];
    y333 = [zz; y333; zz;];
    y444 = [zz; y444; zz;];
    ttt = 0:1/Fs:(0.022 + 1/Fs);
    
    figure;
    hold on;
    plot(ttt,y111);
    plot(ttt,y222);
    plot(ttt,y333);
    plot(ttt,y444);
    allchirp = [y111 y222 y333 y444];
end

      
        
for q=1:10000
        if Mode==0 | Mode==1 | Mode==2
            for n=1:4
                aoSession = daq.createSession('ni');
                aiSession = daq.createSession('ni');
                aiSession.Rate = Fs;
                aiSession.addAnalogInputChannel('Dev1','ai1', 'Voltage'); % scaled output
                aiSession.DurationInSeconds = 0.05;

                if Mode==1
                    aoSession.Rate = 16000;
                    aoSession.addAnalogOutputChannel('Dev1',0:3,'Voltage');
                    aoSession.queueOutputData([wav1 wav2 wav3 wav4]);
                    aoSession.startBackground
                    pause(2);               
                end
                
                if Mode == 0
                    aoSession.Rate = 16000;
                    aoSession.addAnalogOutputChannel('Dev1',0:3,'Voltage');
                    aoSession.queueOutputData([wav1 wav2 wav3 wav4]);
                    aoSession.startBackground
                    pause(7);               
                end

                if Mode==2
                    aoSession.Rate = Fs;
                    aoSession.addAnalogOutputChannel('Dev1',strcat('ao',num2str(n-1)),'Voltage');
                    aoSession.queueOutputData((allchirp(:,n))); 
%                     figure;
%                     plot(tttt,allchirp(:,n));
                    aoSession.startBackground

                    pause(4);
                end                
            end   
        end
        
        if Mode==3
                aoSession = daq.createSession('ni');
                aoSession.Rate = Fs;
                aoSession.addAnalogOutputChannel('Dev1',0:3,'Voltage');              
                aoSession.queueOutputData( gain * [y111 y222 y333 1.5*y444]); 
                aoSession.startBackground
                pause(0.5)
        end      
end
    
 

