function DaqHandler( mode ,loops, Fs, allchirp, pause_time,gain )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
for q = 1:loops
        if mode==0 || mode==1 || mode==2
            for n=1:4
                aoSession = daq.createSession('ni');
                aiSession = daq.createSession('ni');
                aiSession.Rate = Fs;
                aiSession.addAnalogInputChannel('Dev1','ai1', 'Voltage'); % scaled output
                aiSession.DurationInSeconds = 0.05;
                
                if mode==0 || mode==1
                    aoSession.Rate = 16000;
                    aoSession.addAnalogOutputChannel('Dev1',0:3,'Voltage');
                    aoSession.queueOutputData([allchirp(:,1) allchirp(:,2) allchirp(:,3) allchirp(:,4)]);
                    aoSession.startBackground
                    pause(pause_time * 10);               
                else
                    aoSession.Rate = Fs;
                    aoSession.addAnalogOutputChannel('Dev1',strcat('ao',num2str(n-1)),'Voltage');
                    aoSession.queueOutputData((allchirp(:,n))); 
                    aoSession.startBackground
                    pause(pause_time * 4);
                end                
            end   
        end
        
        if mode==3
                aoSession = daq.createSession('ni');
                aoSession.Rate = Fs;
                aoSession.addAnalogOutputChannel('Dev1',0:3,'Voltage');              
                aoSession.queueOutputData( gain * [allchirp(:,1) allchirp(:,2) allchirp(:,3) allchirp(:,4)]); 
                aoSession.startBackground
                pause(pause_time)
        end      
end

end

