function pilots = nonHTPilots(Nsym,z,varargin)
%nonHTPilots Non-HT pilot sequence
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   PILOTS = nonHTPilots(NSYM,Z) returns Non-HT pilots. Null subcarriers
%   are not included therefore PILOTS is sized 4-by-Nsym, where Nsym is the
%   number of symbols within the field.
%
%   NSYM is a scalar specifying the number of symbols within the VHT field.
%
%   Z is a scalar specifying the number of symbols preceding the current
%   field, and is given in the standard as an addition in the pilot
%   polarity sequence subscript, e.g. the 1 in p_{n+1} in IEEE 802.11-2012
%   Eqn 18-22.
%
%   PILOTS = nonHTPilots(NSYM,Z,CHANBW) returns the Non-HT pilots
%   replicated for bandwidths greater than 20 MHz. PILOTS is sized
%   Np-by-Nsym, where Np is the number of occupied pilot subcarriers over
%   the whole channel bandwidth. CHANBW is a string specifying the channel
%   bandwidth.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

if nargin > 2
    chanBW = varargin{1};
else
    chanBW = 'CBW20'; % default
end

n = (0:Nsym-1).'; % Indices of symbols within the field
pilotSeq = nonHTPilotSequence(chanBW);      % IEEE Std 802.11-2012 Eqn 18-24
polaritySeq = pilotPolaritySequence(n+z).'; % IEEE Std 802.11-2012 Eqn 18-25 
pilots = bsxfun(@times,polaritySeq,pilotSeq);

end

function pilotSeq = nonHTPilotSequence(chanBW)

pilots20Mhz = [1 1 1 -1].'; % IEEE Std 802.11-2012 Eqn 18-24.

if (strcmp(chanBW,'CBW5') || strcmp(chanBW,'CBW10') || strcmp(chanBW,'CBW20')) 
    % Same FFT length for 5/10/20 MHz
    num20 = 1;
else
    num20 = real(str2double(chanBW(4:end)))/20;
end

% Replicate pilots over 20 MHz segments
pilotSeq = repmat(pilots20Mhz,num20,1);

end