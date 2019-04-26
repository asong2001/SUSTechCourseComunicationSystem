function Q = helperSpatialExpansionMatrix(cfg)
% helperSpatialExpansionMatrix Return a spatial expansion matrix
%
%   Q = helperSpatialExpansionMatrix(CFGFORMAT) returns a spatial expansion
%   matrix for the specified format configuration object, CFGFORMAT.
%
%   Q is a complex matrix sized Nst-by-Nsts-by-Nt. Nst is the number of
%   occupied subcarriers, Nsts is the number of space-time streams, and Nt
%   is the number of transmit antennas.
%
%   CFGFORMAT is the format configuration object of type <a
%   href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> or <a
%   href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, which specifies
%   the parameters for the VHT or HT-Mixed formats respectively.
%
%   Example: Spatial expansion for a VHT format configuration. 
%
%   cfgVHT = wlanVHTConfig;
%   cfgVHT.SpatialMapping = 'Custom';
%   Q = helperSpatialExpansionMatrix(cfgVHT);
%   cfgVHT.SpatialMappingMatrix = Q;
%
%   See also wlanVHTConfig, wlanHTConfig.

%   Copyright 2015 The MathWorks, Inc.

validateattributes(cfg,{'wlanVHTConfig','wlanHTConfig'}, ...
    {'scalar'}, mfilename,'format configuration object');

NumSTS = cfg.NumSpaceTimeStreams;
NumTx = cfg.NumTransmitAntennas;

csd = helperCyclicShiftDelay(NumTx,helperSampleRate(cfg),'VHT');
Nfft = helperFFTLength(cfg);
[dataInd,pilotInd] = helperSubcarrierIndices(cfg,'VHT');

% Calculate MCSD matrix (Std 802.11-2012 Section 20.3.11.2)
n = (1:Nfft)-Nfft/2-1;
phaseShift = exp(-1i*2*pi*csd*n/Nfft).';
Mcsd = permute(phaseShift,[1 3 2]);

% Calculate D matrix (Std 802.11-2012 Section 20.3.11.2)
base = eye(NumSTS);
D = zeros(NumTx,NumSTS);
for d=1:NumTx
   D(d,:) = base(mod(d-1,NumSTS)+1,:); 
end
D = sqrt(NumSTS/NumTx)*D;

% Calculate spatial expansion matrix for occupied subcarriers
fullMat = bsxfun(@times,Mcsd,permute(D,[3 2 1]));
Q = fullMat(sort([dataInd; pilotInd]),:,:);
end