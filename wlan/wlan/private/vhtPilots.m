function pilots = vhtPilots(Nsym,z,chanBW,numSTS)
%vhtPilots VHT pilot sequence
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   PILOTS = vhtPilots(NSYM,Z,CHANBW,NUMSTS) returns VHT pilots.
%
%   PILOTS is Nsp-by-Nsym-by-Nsts, where Nsp is the number of pilot
%   subcarriers, Nsym is the number of OFDM symbols, and Nsts is the number
%   of space time streams.
%
%   NSYM is a scalar specifying the number of symbols within the VHT field.
%
%   Z is a scalar specifying the number of symbols preceding the current
%   field, and is given in the standard as an addition in the pilot
%   polarity sequence subscript, e.g. the 4 in p_{n+4} in IEEE
%   802.11ac-2013 Eqn 22-95.
%
%   CHANBW is the channel bandwidth string.
%
%   NUMSTS is the number of space-time streams.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

n = (0:Nsym-1).'; % Indices of symbols within the field
pilotSeq = vhtPilotSequence(chanBW,numSTS,n); % IEEE Std 802.11ac-2013 Section 22.3.10.10
polaritySeq = pilotPolaritySequence(n+z).';   % IEEE Std 802.11-2012 Eqn 18-25 
pilots = bsxfun(@times,polaritySeq,pilotSeq); 
end

% VHT Pilot sequence defined in IEEE Std 802.11ac-2013 Section 22.3.10.10.
% Returns a matrix sized Nsp-by-NSym-by-Nsts, where Nsp is the number of
% pilot subcarriers, NSym is the number of OFDM symbols, and NSTS is the
% number of space time streams.
function pilotSeq = vhtPilotSequence(chanBW,numSTS,n)
% Force n to be a column array for codegen
if ~iscolumn(n)
    nC = n.';
else
    nC = n;
end

Psi80MHz = [1; 1; 1; -1; -1; 1; 1; 1]; % IEEE Std 802.11ac-2013 Table 22-21
switch chanBW
    case {'CBW20','CBW40'}
        % Uses HT pilot sequence with a single space-time stream
        singleSTSPilots = htPilotSequence(chanBW,1,nC);
    case {'CBW80','CBW80+80'}
        P = numel(Psi80MHz);
        idx = mod(bsxfun(@plus,nC(:,1),(0:P-1)),P)+1;
        singleSTSPilots = Psi80MHz(idx.');
    otherwise % {'CBW160'}
        Psi160MHz = [Psi80MHz; Psi80MHz];
        P = numel(Psi160MHz);
        idx = mod(bsxfun(@plus,nC(:,1),(0:P-1)),P)+1;
        singleSTSPilots = Psi160MHz(idx.');
end
pilotSeq = repmat(singleSTSPilots,1,1,numSTS); % Same pilots per STS
end