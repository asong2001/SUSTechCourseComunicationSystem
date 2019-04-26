function noiseVar = vhtBeamformingNoiseVariance(noisePower,cfgVHT)
% vhtBeamformingNoiseVariance returns the noise variance after OFDM
% demodulation for a given noise power.
%
%   NOISEVAR = vhtBeamformingNoiseVariance(NOISEPOWER,CFGVHT) returns the
%   nosie variance after OFDM demodulation for a given noise power,
%   NOISEPOWER, and VHT format configuration object, CFGFORMAT.
%
%   NOISEPOWER is the noise power applied to each receive antenna in dBW.
%
%   CFGVHT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>.

% Copyright 2015 The MathWorks, Inc.

%#codegen

% An estimate of the noise power after OFDM demodulation is required to
% perform MMSE equalization on the received OFDM symbols. The noise
% variance after demodulation is calculated and is used during field
% recovery. The noise variance after OFDM demodulation in VHT fields is a
% scaled version of the applied noise power. This scaling accounts for:
%
% * The noise energy in unused subcarriers which are removed during
% demodulation
% * Scaling by the number of space-time streams during demodulation 

% Get the number of occupied subcarriers in VHT fields
[data,pilots] = helperSubcarrierIndices(cfgVHT,'VHT');
Nsd = numel(data);   % Number of data subcarriers
Nsp = numel(pilots); % Number of pilot subcarriers
Nst = Nsd+Nsp;
Nfft = helperFFTLength(cfgVHT); % FFT length

% Calculate VHT noise variance
noiseVar = 10^(noisePower/10)*cfgVHT.NumSpaceTimeStreams*(Nst/Nfft);

end