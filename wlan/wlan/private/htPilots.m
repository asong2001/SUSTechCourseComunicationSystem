function pilots = htPilots(Nsym,z,chanBW,Nsts)
%htPilots HT pilot sequence
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   PILOTS = htPilots(NSYM,Z,CHANBW,NSTS) returns HT pilots.
%
%   PILOTS is Nsp-by-NSym-by-NSTS, where Nsp is the number of pilot
%   subcarriers, NSym is the number of OFDM symbols, and NSTS is the number
%   of space time streams.
%
%   NSYM is a scalar specifying the number of symbols within the HT field.
%
%   Z is a scalar specifying the number of symbols preceding the current
%   field, and is given in the standard as an addition in the pilot
%   polarity sequence subscript, e.g. the 1 in p_{n+1} in IEEE 802.11-2012
%   Eqn 20-17.
%
%   CHANBW is the channel bandwidth string.
%
%   NSTS is the number of space-time streams.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

n = (0:Nsym-1).'; % Indices of symbols within the field
pilotSeq = htPilotSequence(chanBW,Nsts,n);  % IEEE Std 802.11-2012 Section 20.3.11.10
polaritySeq = pilotPolaritySequence(n+z).'; % IEEE Std 802.11-2012 Eqn 18-25 
pilots = bsxfun(@times,polaritySeq,pilotSeq);
end
