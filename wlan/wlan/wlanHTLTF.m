function y = wlanHTLTF(cfgHT)
%wlanHTLTF HT Long Training Field (HT-LTF)
%
%   Y = wlanHTLTF(CFGHT) generates the HT Long Training Field (HT-LTF)
%   time-domain waveform for the HT-Mixed transmission format.
%
%   Y is the time-domain HT-LTF signal. It is a complex matrix of size
%   Ns-by-Nt where Ns represents the number of time-domain samples and
%   Nt represents the number of transmit antennas.
%
%   CFGHT is the format configuration object of type <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a> which
%   specifies the parameters for the HT-Mixed format.
%
%   Example: 
%   %  Generate the HT-LTF waveform for a HT 40MHz transmission format
%
%      cfgHT = wlanHTConfig;                % Format configuration
%      cfgHT.ChannelBandwidth = 'CBW40';    % Set to 40MHz
%      hltfOut = wlanHTLTF(cfgHT);
%
%   See also wlanHTConfig, wlanLLTF, wlanHTData, wlanHTLTFDemodulate.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

validateattributes(cfgHT, {'wlanHTConfig'}, {'scalar'}, mfilename, ...
                   'HT-Mixed format configuration object');
validateConfig(cfgHT, 'SMapping'); 

chanBW = cfgHT.ChannelBandwidth;
numSTS = cfgHT.NumSpaceTimeStreams;
numTx = cfgHT.NumTransmitAntennas;
spatialMapMtx = cfgHT.SpatialMappingMatrix;
if inESSMode(cfgHT)
    numESS = cfgHT.NumExtensionStreams;
else
    numESS = 0;
end

% Get OFDM parameters
cfgOFDM = wlanGetOFDMConfig(chanBW, 'Long', 'HT', numSTS);
FFTLen      = cfgOFDM.FFTLength;
CPLen       = cfgOFDM.CyclicPrefixLength;
dataIdx     = cfgOFDM.DataIndices;
pilotIdx    = cfgOFDM.PilotIndices;
chanBWInMHz = FFTLen/64 * 20;

% HT training fields are subset of VHT
[HTLTF, Phtltf, Ndltf, Neltf] = vhtltfSequence(chanBW, numSTS, numESS);
Nltf = Ndltf + Neltf;

htltfToneRotated = HTLTF .* cfgOFDM.CarrierRotations;

% Define HTLTF and output variable sizes
htltfLen = FFTLen + CPLen;
y = complex(zeros(htltfLen*Nltf, numTx));

% Generate HT-Data LTFs
htltfSTS = complex(zeros(FFTLen, numSTS));
Pd = Phtltf(1:numSTS,1:Ndltf);
csh = getCyclicShiftVal('VHT', numSTS, chanBWInMHz); 
for i = 1:Ndltf
    
    htltfSTS(dataIdx,:) = repmat(htltfToneRotated(dataIdx), 1, numSTS) .* ...
        repmat(Pd(:,i).', length(dataIdx), 1);
    
    htltfSTS(pilotIdx,:)= repmat(htltfToneRotated(pilotIdx), 1, numSTS) .* ...
        repmat(Pd(:,i).', length(pilotIdx), 1);
    
    % Cyclic shift addition
    % The cyclic shift is applied per stream.
    htltfCycShift = wlanCyclicShift(htltfSTS, csh, FFTLen, 'Tx');
        
    % Spatial mapping
    if strcmp(cfgHT.SpatialMapping, 'Custom')
        if ismatrix(spatialMapMtx)
            if isscalar(spatialMapMtx) || isvector(spatialMapMtx)
                htlfSpatialMapped = wlanSpatialMapping(htltfCycShift, ...
                    cfgHT.SpatialMapping, numTx, spatialMapMtx);
            else
                htlfSpatialMapped = wlanSpatialMapping(htltfCycShift, ...
                    cfgHT.SpatialMapping, numTx, spatialMapMtx(1:numSTS, :));
            end
        else % 3D
            htlfSpatialMapped = wlanSpatialMapping( ...
                htltfCycShift, cfgHT.SpatialMapping, ...
                numTx, spatialMapMtx(:, 1:numSTS, :));
        end
    else
        htlfSpatialMapped = wlanSpatialMapping(htltfCycShift, ...
            cfgHT.SpatialMapping, numTx, spatialMapMtx);
    end
    
    % OFDM modulation
    modOut = ifft(ifftshift(htlfSpatialMapped, 1), [], 1);
    tmp = [modOut(end - CPLen + 1:end,:); modOut];

    y((i-1)*htltfLen + (1:htltfLen).',:) = tmp * cfgOFDM.NormalizationFactor; 
end

% Append the HT-Extension LTFs as well, if specified
if numESS>0
    htltfESS = complex(zeros(FFTLen, numESS));
    Pe = Phtltf(1:numESS,1:Neltf);
    csh = getCyclicShiftVal('VHT', numESS, chanBWInMHz);

    for i = 1:Neltf
        
        htltfESS(dataIdx,:) = repmat(htltfToneRotated(dataIdx),1, ...
            numESS).*repmat(Pe(:,i).', length(dataIdx), 1);
        
        htltfESS(pilotIdx,:)= repmat(htltfToneRotated(pilotIdx),1, ...
            numESS).*repmat(Pe(:,i).', length(pilotIdx), 1);
        
        % Cyclic shift addition
        % The cyclic shift is applied per stream.
        htltfCycShift = wlanCyclicShift(htltfESS, csh, FFTLen, 'Tx');
        
        % Spatial mapping
        if ismatrix(spatialMapMtx)
            if isscalar(spatialMapMtx) || isvector(spatialMapMtx)
                htlfSpatialMapped = wlanSpatialMapping(htltfCycShift, ...
                    cfgHT.SpatialMapping, numTx, spatialMapMtx);
            else
                htlfSpatialMapped = wlanSpatialMapping(htltfCycShift, ...
                    cfgHT.SpatialMapping, numTx, spatialMapMtx(numSTS+1:end, :) );
            end
        else % 3D
            htlfSpatialMapped = wlanSpatialMapping(htltfCycShift, ...
                cfgHT.SpatialMapping, numTx, spatialMapMtx(:, numSTS+1:end, :) );
        end
        
        % OFDM modulation
        modOut = ifft(ifftshift(htlfSpatialMapped, 1), [], 1);
        tmp = [modOut(end - CPLen + 1:end,:); modOut];

        y(htltfLen*Ndltf + (i-1)*htltfLen + (1:htltfLen).',:) = ...
            tmp * cfgOFDM.NormalizationFactor * sqrt(numSTS/numESS);
    end    
end

end

% [EOF]
