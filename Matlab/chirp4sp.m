function [ sig ] = chirp4sp( t,start,time,stop )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    sig = chirp(t,start,time,stop);
    sig=sig';
end

