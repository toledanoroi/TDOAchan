% ni-daq temprature couple measurement and averaging
% Authors: Roi Toledano and Yarden Avraham, 

% just for debug 
devices = daq.getDevices;
my_dev = devices(1);

% create Temprature measurements session
s = daq.createSession('ni');
s.Rate = 10;
s.addAnalogInputChannel('Dev1',0, 'Thermocouple');

% define channel parameters:
%            Units: [ Celsius | Fahrenheit | Kelvin | Rankine ]
% ThermocoupleType: [ Unknown | J | K | N | R | S | T | B | E ]
%         Coupling: [ DC | AC ]
%   TerminalConfig: [ Differential | SingleEnded | SingleEndedNonReferenced | PseudoDifferential ]
%            Range: 0 to +750 Celsius
%             Name: {}
tc = s.Channels(1);
tc.Units = 'Celsius';
% tc.ThermocoupleType = 'K';
[data,time] = s.startForeground();

wait(s,30);
stop(s);

plot(time, data)
xlabel('Time (secs)');
ylabel('Temperature (Kelvin)');

% calculate statistical parameters
mu = mean(data);
variance = var(data);
stdeviation = std(data);


save('Temprature.mat','mu','variance','stdeviation');

