function [ sig ] = UpAndDownchirp4sp( t,low,time,high,mode)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here  
    if mode == 'up'
        y1=chirp(t,low,time,high);
        y11=chirp(t,high,time,low,'linear', acos(y1(length(y1)-1)) * 180/pi );
        y111=[y1,y11];
    elseif mode == 'down'
        y1 = chirp(t,high,time,low);
        y11=chirp(t,low,time,high,'linear', acos(y1(length(y1)-1)) * 180/pi );
        y111=[y1,y11];
    end
    sig=y111';
end

