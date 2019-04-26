function rxWaveform = vhtSigRecGenerateWaveform(cfgVHTTx,numRx, ...
    delayProfile,noisePower,cfo,numTxPkt,idleTime)
%vhtSigRecGenerateWaveform Featured example helper function
%
%   Generates an impaired multi-packet VHT waveform

%   Copyright 2015 The MathWorks, Inc.

rngState = rng(120); % Set random state
txPSDUPerPkt = randi([0,1],8*cfgVHTTx.PSDULength,1,'int8');

% Generate multi-packet waveform with PSDU repeatedly applied. With each
% packet containing the same PSDU, this enables us to perform packet error
% checking after recovering the packets.
txWave = wlanWaveformGenerator(txPSDUPerPkt,cfgVHTTx, ...
    'NumPackets',numTxPkt,'IdleTime',idleTime);

% Pass through a TGac fading channel with leading zeros
sr = helperSampleRate(cfgVHTTx.ChannelBandwidth);
TGacChan = wlanTGacChannel( ...
    'SampleRate',             sr, ...
    'DelayProfile',           delayProfile, ...             
    'ChannelBandwidth',       cfgVHTTx.ChannelBandwidth, ...
    'NumTransmitAntennas',    cfgVHTTx.NumTransmitAntennas, ...
    'NumReceiveAntennas',     numRx, ...
    'LargeScaleFadingEffect', 'None'); % No path loss or shadowing
chanOut = step(TGacChan, ...
    [zeros(round(idleTime*sr),cfgVHTTx.NumTransmitAntennas); txWave]);

% Add carrier frequency offset
chanOut = helperFrequencyOffset(chanOut, sr, cfo);

% Add noise
AWGN = comm.AWGNChannel( ...
    'NoiseMethod', 'Variance', ...
    'Variance',    10^(noisePower/10));
rxWaveform = step(AWGN,chanOut);
rng(rngState); % Restore random state

end

