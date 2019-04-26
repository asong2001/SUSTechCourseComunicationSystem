function y = wlanVHTLTF(cfgVHT)
%WLANVHTLTF VHT Long Training Field (VHT-LTF)
% 
%   Y = wlanVHTLTF(CFGVHT) generates the VHT Long Training Field (VHT-LTF)
%   time-domain signal for the VHT transmission format.
%
%   Y is the time-domain VHT-LTF signal. It is a complex matrix of size
%   Ns-by-Nt where Ns represents the number of time-domain samples and Nt
%   represents the number of transmit antennas.
%
%   CFGVHT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> which
%   specifies the parameters for the VHT format.
% 
%   Example: 
%   %  Generate the VHT-LTF signal for a VHT 80MHz transmission format
% 
%      cfgVHT = wlanVHTConfig;                % Format configuration
%      cfgVHT.ChannelBandwidth = 'CBW80';     % Set to 80MHz
%      vhtltfOut = wlanVHTLTF(cfgVHT);
% 
%   See also wlanVHTConfig, wlanLLTF, wlanVHTSTF, wlanVHTData,
%   wlanVHTLTFDemodulate.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, mfilename, ...
                   'VHT format configuration object');
validateConfig(cfgVHT, 'SMapping');

chanBW = cfgVHT.ChannelBandwidth;
numSTSTotal = sum(cfgVHT.NumSpaceTimeStreams); 
numTx = cfgVHT.NumTransmitAntennas; 
   
% Get OFDM parameters
cfgOFDM = wlanGetOFDMConfig(chanBW, 'Long', 'VHT', numSTSTotal);
FFTLen = cfgOFDM.FFTLength;
CPLen  = cfgOFDM.CyclicPrefixLength;

% Get VHT-LTF sequences
[VLTF, Pvhtltf, Nltf] = vhtltfSequence(chanBW, numSTSTotal);

vltfToneRotated = VLTF .* cfgOFDM.CarrierRotations;

% Define VLTF and output variable sizes
vhtltfSTS = complex(zeros(FFTLen, numSTSTotal));
vltfLength = FFTLen + CPLen;
y = complex(zeros(vltfLength*Nltf, numTx));

% P and R matrices as per IEEE Std 802.11ac-2013 Sec. 22.3.8.3.5
P = Pvhtltf(1:numSTSTotal,1:Nltf);
R = Pvhtltf(1,1:Nltf);

csh = getCyclicShiftVal('VHT', numSTSTotal, FFTLen/64*20);

% Generate and modulate each VHT-LTF symbol
for i = 1:Nltf
    
    vhtltfSTS(cfgOFDM.DataIndices,:) = bsxfun(@times, ...
        vltfToneRotated(cfgOFDM.DataIndices), P(:, i).');    
    vhtltfSTS(cfgOFDM.PilotIndices,:)= ...
        repmat(vltfToneRotated(cfgOFDM.PilotIndices), 1, numSTSTotal) * R(1, i);
    
    % Cyclic shift addition.
    % The cyclic shift is applied per user per stream.
    vltfCycShift =  wlanCyclicShift(vhtltfSTS, csh, FFTLen, 'Tx');
    
    % Spatial mapping
    vltfSpatialMapped = wlanSpatialMapping(vltfCycShift, cfgVHT.SpatialMapping, ...
        cfgVHT.NumTransmitAntennas, cfgVHT.SpatialMappingMatrix);

    % OFDM modulation
    modOut = ifft(ifftshift(vltfSpatialMapped, 1), [], 1); 
    tmp = [modOut(end - CPLen + 1:end,:); modOut];
    y((i-1)*vltfLength + (1:vltfLength),:) = ...
                                tmp * cfgOFDM.NormalizationFactor;
                            
end

end

% [EOF]
