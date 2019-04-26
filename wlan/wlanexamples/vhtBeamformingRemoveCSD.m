function y = vhtBeamformingRemoveCSD(x,cfgFormat)
% vhtBeamformingRemoveCSD Remove effect of transmitter cyclic shifts
%
%   Y = removeChanEstCyclicShifts(X,CFGFORMAT) returns a channel estimate
%   matrix without the effect of cyclic shifts applied at the transmitter.
%
%   Y is a complex Nst-by-Nsts-by-Nr array containing the estimated channel
%   at data and pilot subcarriers without the effects of transmitter cyclic
%   shifts. Nst is the number of occupied subcarriers, Nsts is the number
%   of space-time streams and Nr is the number of receive antennas.
%
%   X is a complex Nst-by-Nsts-by-Nr array containing the estimated channel
%   at data and pilot subcarriers with the effects of transmitter cyclic
%   shifts.
%
%   CFGFORMAT is the format configuration object of type <a
%   href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> or <a
%   href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, which specifies
%   the parameters for the VHT or HT-Mixed formats respectively.

%   Copyright 2015 The MathWorks, Inc.

csd = helperCyclicShiftDelay(cfgFormat.NumSpaceTimeStreams, ...
    helperSampleRate(cfgFormat),'VHT');
Nfft = helperFFTLength(cfgFormat);
[dataIndices,pilotIndices] = helperSubcarrierIndices(cfgFormat,'VHT');
occupiedIndices = sort([dataIndices;pilotIndices]);

% Remove cyclic shift (note the negative)
k = (1:Nfft)-Nfft/2-1;
phaseShift = exp(-1i*2*pi*(-csd)*k/Nfft).';
y = bsxfun(@times,x,phaseShift(occupiedIndices,:));    
end