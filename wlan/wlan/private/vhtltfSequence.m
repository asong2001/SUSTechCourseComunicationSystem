function [ltfSC,varargout] = vhtltfSequence(chanBW,numSTS,varargin)
%vhtltfSequence HT-LTF, VHT-LTF subcarrier sequence and related parameters
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   [LTFSYM,P,NUMDLTF] = vhtltfSequence(CHANBW,NUMSTS) returns the
%   subcarrier values for HT-LTF/VHT-LTF, P matrix and the number of LTF
%   symbols, for a given channel bandwidth CHANBW and number of space-time
%   streams NUMSTS, for HT-Mixed and VHT formats.
%
%   [LTFSYM,P,NUMDLTF,NUMELTF] = vhtltfSequence(CHANBW,NUMSTS,NUMESS)
%   returns the subcarrier values for the HT-LTF, P matrix and number of
%   HT-LTF symbols (both data and extension) for a given channel bandwidth
%   CHANBW, number of space-time streams NUMSTS and number of extension
%   streams NUMESS, for the HT-Mixed format.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(2,3);
nargoutchk(0,4);
if nargin>2
    htMode = true;
    numESS = varargin{1}; % in range [0, 3]
else 
    htMode = false;
    numESS = 0;
end

% Subcarrier values for the VHT-LTF symbol including guard bands and DC
[ltfLeft, ltfRight] = lltfSequence(); % Sequences based on lltf

switch chanBW
    case {'CBW20'}
        % Same for HT and VHT
        ltfSC = [zeros(4,1); 1; 1; ltfLeft; 0; ltfRight;-1;-1; zeros(3,1)];
    case {'CBW40'}
        % Same for HT and VHT
        ltfSC = [ zeros(6,1); ltfLeft; 1; ltfRight;-1;-1;-1; 1; 0; 0; 0; ...
            -1;1; 1;-1; ltfLeft; 1; ltfRight;  zeros(5,1)];
    case {'CBW80'}
        % VHT only
        ltfSC = [zeros(6,1); ltfLeft; 1; ltfRight;-1;-1;-1; 1; 1;-1; 1; ...
            -1; 1; 1;-1; ltfLeft; 1; ltfRight; 1;-1; 1;-1; 0; 0; 0; ...
            1;-1;-1; 1; ltfLeft; 1; ltfRight;-1;-1;-1; 1; 1;-1; 1; ...
            -1; 1; 1;-1; ltfLeft; 1; ltfRight; zeros(5,1)];
    otherwise % {'CBW160'}
        % VHT only
        vltf80 = [ltfLeft; 1; ltfRight;-1;-1;-1; 1; 1;-1; 1;-1; 1; 1;-1; ...
            ltfLeft; 1; ltfRight; 1;-1; 1;-1; 0; 0; 0; 1;-1;-1; 1; ...
            ltfLeft; 1; ltfRight;-1;-1;-1; 1; 1;-1; 1;-1; 1; 1;-1; ...
            ltfLeft; 1; ltfRight];
        ltfSC = [zeros(6,1); vltf80; 0; 0; 0; 0; 0; 0; 0; ...
            0; 0; 0; 0; vltf80; zeros(5,1)];
end

if (nargout>1)
    % Mapping matrix
    % IEEE Std 802.11-2012 Eqn 20-27
    Phtltf = [1 -1 1 1; 1 1 -1 1; 1 1 1 -1; -1 1 1 1];
    if (sum(numSTS) <= 4) % 1-4
        P = Phtltf;
    elseif sum(numSTS) > 4 && sum(numSTS) <= 6  % numSTS = {5,6}
        w = exp(-2*pi*1i/6);
        % IEEE Std 802.11ac-2013 Eqn 22-44
        P = [1 -1   1    1    1    -1; ...
            1 -w   w^2  w^3  w^4  -w^5; ...
            1 -w^2 w^4  w^6  w^8  -w^10; ...
            1 -w^3 w^6  w^9  w^12 -w^15; ...
            1 -w^4 w^8  w^12 w^16 -w^20; ...
            1 -w^5 w^10 w^15 w^20 -w^25 ];
    else    % numSTS ={7,8}
        P = [Phtltf Phtltf; Phtltf -Phtltf];
    end
    varargout{1} = P;
end

if (nargout>2)
    %Number of HT/VHT LTF symbols
    if htMode % HT-mixed
        % Includes ELTFs in addition to DLTFs
        NhtdltfTable = [1 2 4 4]; % for only HTDLTFs
        NhteltfTable = [0 1 2 4]; % for only HTELTFS
        numDLTFSym = NhtdltfTable(numSTS);
        numELTFSym = NhteltfTable(numESS+1);
        if (nargout==4)
            varargout{3} = numELTFSym;
        end
    else % VHT
        % Number of VHT-LTF symbols required for different numbers of
        % space time streams: IEEE Std 802.11ac-2013 Table 22-13, page 264.
        NvhtltfTable = [1 2 4 4 6 6 8 8];
        numDLTFSym = NvhtltfTable(numSTS); % Number of VHTLTF symbols
    end
    varargout{2} = numDLTFSym;
end

end