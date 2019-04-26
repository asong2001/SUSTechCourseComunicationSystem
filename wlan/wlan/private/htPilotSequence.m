function pilots = htPilotSequence(chanBW,NSTS,n)
%htPilotSequence HT pilot sequence
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   PILOTS = htPilotSequence(CHANBW,NSTS,N) returns the HT pilot sequence
%   specified in Std 802.11-2012 Section 20.3.11.10. PILOTS is sized
%   Nsp-by-Nsym-by-Nsts, where Nsp is the number of pilots, Nsym is the
%   number of OFDM symbols, and Nsts is the number of space time streams.
%
%   CHANBW is the channel bandwidth string.
%   NSTS is the total number of space time streams (1:4)
%   N is the index of the OFDM symbol (0:...)

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Force n to be a column array for codegen
if ~iscolumn(n)
    nC = n.';
else
    nC = n;
end

% Tone packing for different CBW for a single stream, single symbol Pilot
% packing
switch chanBW
    case {'CBW20'}
        % Per row for one mode of NSTS and iSTS - Table 20-19 Each row
        % holds the pilots for all subcarriers per symbol
        pBasePilotMatrix = [ ...
            1 1 1 -1 ; ... % NSTS:1 iSTS:1
            1 1 -1 -1; ... % NSTS:2 iSTS:1
            1 -1 -1 1; ... % NSTS:2 iSTS:2
            1 1 -1 -1; ... % NSTS:3 iSTS:1
            1 -1 1 -1;...  % NSTS:3 iSTS:2
            -1 1 1 -1; ... % NSTS:3 iSTS:3
            1 1 1 -1 ; ... % NSTS:4 iSTS:1
            1 1 -1 1 ; ... % NSTS:4 iSTS:2
            1 -1 1 1 ; ... % NSTS:4 iSTS:3
            -1 1 1 1 ];    % NSTS:4 iSTS:4
        pPilotMod = 4;               % Eqn. 20-54, pg 1720
    otherwise % case {'CBW40'}
        % Per row for one mode of NSTS and iSTS - Table 20-20
        pBasePilotMatrix = [ ...
            1 1 1 -1 -1 1  ; ... % NSTS:1 iSTS:1
            1 1 -1 -1 -1 -1; ... % NSTS:2 iSTS:1
            1 1 1 -1 1 1   ; ... % NSTS:2 iSTS:2
            1 1 -1 -1 -1 -1; ... % NSTS:3 iSTS:1
            1 1 1 -1 1 1   ; ... % NSTS:3 iSTS:2
            1 -1 1 -1 -1 1 ; ... % NSTS:3 iSTS:3
            1 1 -1 -1 -1 -1; ... % NSTS:4 iSTS:1
            1 1 1 -1 1 1   ; ... % NSTS:4 iSTS:2
            1 -1 1 -1 -1 1 ; ... % NSTS:4 iSTS:3
            -1 1 1 1 -1 1];      % NSTS:4 iSTS:4
        pPilotMod = 6;              % Eqn. 20-55, pg 1720
end
iSTS = (1:NSTS).'; % Indices of space-time streams 

stIdxVec = [1 2 4 7]; % Row offset for each NSTS
idx = (mod(bsxfun(@plus,nC(:,1),0:pPilotMod-1),pPilotMod)+1).';
pilots = pBasePilotMatrix(stIdxVec(NSTS)+iSTS-1,idx(:));

% Reshape from STS-by-K-by-N to K-by-N-by-STS
pilots = permute(reshape(pilots,numel(iSTS),pPilotMod,numel(nC)),[2 3 1]);
end
