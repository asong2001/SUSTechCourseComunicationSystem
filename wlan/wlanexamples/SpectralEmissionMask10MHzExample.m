%% 802.11p(TM) Spectral Emission Mask Testing
%
% This example shows how to perform spectrum emission mask tests for an
% IEEE(R) 802.11p(TM) transmitted waveform.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% IEEE 802.11p [ <#11 2> ] is an approved amendment to the IEEE 802.11
% standard to enable support for wireless access in vehicular environments
% (WAVE). Using the half-clocked mode with a 10 MHz channel bandwidth, it
% operates at the 5.85-5.925 GHz bands for which additional spectral
% emission masks are defined [ Annex D of <#11 1> ].
%
% This example shows how spectral mask measurements can be performed on a
% transmitted waveform. The waveform is generated with WLAN System
% Toolbox(TM) for simplicity, but a waveform captured with a spectrum
% analyzer could be used as well.
%
% A waveform consisting of three 10 MHz 802.11p packets separated by a 32
% microsecond gap is generated. Random data is used for each packet and
% 16QAM modulation is used. The baseband waveform is upsampled and filtered
% to reduce the out of band emissions to meet the spectral mask
% requirement. A high power amplifier (HPA) model is used, which introduces
% inband distortion and spectral regrowth. The spectral emission mask
% measurement is performed on the upsampled waveform after the high power
% amplifier modeling. The test schematic is illustrated in the following
% diagram:
%
% <<SpectralEmissionMaskDiagram.png>>

%% IEEE 802.11p non-HT Packet Configuration
% In this example, an IEEE 802.11p waveform consisting of multiple non-HT
% format packets is generated. The format specific configuration of the
% non-HT waveform is described using a non-HT format configuration object.
% The object is created using the <matlab:doc('wlanNonHTConfig')
% wlanNonHTConfig> function. In this example, the object is configured for
% a 10 MHz bandwidth operation as used by IEEE 802.11p.

cfgNHT = wlanNonHTConfig;          % Create packet configuration
cfgNHT.ChannelBandwidth = 'CBW10'; % 10 MHz
cfgNHT.MCS = 4;                    % Modulation 16QAM, rate-1/2
cfgNHT.PSDULength = 1000;          % PSDU length in bytes

%% Baseband Waveform Generation
% The waveform generator can be configured to generate one or more packets
% and add an idle time between each packet. In this example three packets
% with a 32 microsecond idle period will be created. Random bits for all
% packets |data| are created and passed as an argument to
% <matlab:doc('wlanWaveformGenerator') wlanWaveformGenerator> along with
% the non-HT packet configuration object |cfgNHT| and additional waveform
% generation parameters. |cfgNHT| configures the waveform generator to
% create the 802.11p Non-HT waveform.

% Set random stream for repeatability of results
s = rng(98765);

% Generate a multi-packet waveform
idleTime   = 32e-6;     % 32 microsecond idle time between packets
numPackets = 3;         % Generate 3 packets

% Create random data; PSDULength is in bytes
data = randi([0 1], cfgNHT.PSDULength*8*numPackets, 1);

genWaveform = wlanWaveformGenerator(data, cfgNHT, ...
                'NumPackets', numPackets,...
                'IdleTime', idleTime);

% Get the sampling rate of the waveform
fs = helperSampleRate(cfgNHT);
disp(['Baseband sampling rate: ' num2str(fs/1e6) ' Msps']);

%% Oversampling and Filtering
% Spectral filtering is used to reduce the out of band spectral emissions
% owing to the implicit rectangular pulse shaping in the OFDM modulation,
% and spectral regrowth caused by the high power amplifier in an RF chain.
% To model the effect of a high power amplifier on the waveform and view
% the out of band spectral emissions the waveform must be oversampled.
% Oversampling requires an interpolation filter to remove spectral images
% caused by upsampling. In this example the waveform is oversampled with an
% interpolation filter which also acts as a spectral filter. This allows
% the waveform to meet spectral mask requirements. The waveform is
% oversampled and filtered using <matlab:doc('dsp.FIRInterpolator')
% dsp.FIRInterpolator>.

% Oversample the waveform
osf = 3;         % Oversampling factor
filterLen = 100; % Filter length
r = 50;          % Design parameter for Chebyshev window (attenuation, dB)

% Generate filter coefficients and interpolate
coeffs = osf.*firnyquist(filterLen, osf, chebwin(filterLen+1, r)); 
coeffs = coeffs(1:end-1);   % Remove trailing zero
FIRINTERP = dsp.FIRInterpolator(osf, 'Numerator', coeffs);
filtWaveform = step(FIRINTERP, [genWaveform; zeros(filterLen/2,1)]);

% Plot the magnitude and phase response of the filter applied after
% oversampling
h = fvtool(FIRINTERP);
h.Analysis = 'freq';           % Plot magnitude and phase responses
h.FS = osf*fs;                 % Set sampling rate
h.NormalizedFrequency = 'off'; % Plot responses against frequency

%% High Power Amplifier Modeling
% Within an RF chain, the high power amplifier is a necessary component
% that also introduces nonlinear behavior in the form of inband distortion
% and spectral regrowth. The Rapp model is used to simulate power
% amplifiers in 802.11. The Rapp model causes AM/AM distortion and is
% modeled with <matlab:doc('comm.MemorylessNonlinearity')
% comm.MemorylessNonlinearity>. The high power amplifier is backed-off to
% operate below the saturation point to reduce distortion. The backoff is
% controlled by the variable |hpaBackoff|.

hpaBackoff = 6; % dB

% Create and configure a memoryless nonlinearity to model the amplifier
HPA = comm.MemorylessNonlinearity;
HPA.Method = 'Rapp model';
HPA.Smoothness = 3;             % p parameter
HPA.LinearGain = -hpaBackoff;   % dB

% Apply the model to the transmit waveform
txWaveform = step(HPA, filtWaveform);

%% Transmit Spectrum Emission Mask Measurement
%
% Stations are classified according to the allowed maximum transmit powers
% (in mW). For the four different classes of stations, four different
% spectral emission masks are defined [ Annex D of <#11 1>]. The spectral
% masks are defined relative to the peak power spectral density (PSD).
%
% In this example the spectrum emission mask of the transmitted
% waveform after high power amplifier modeling is measured for a Class
% A station.

% IEEE Std 802.11-2012 Annex D.2.3, Table D-5: Class A STA
dBrLimits = [-40  -40 -28 -20  -10 0   0  -10 -20 -28 -40 -40];  
fLimits   = [-Inf -15 -10 -5.5 -5 -4.5 4.5 5  5.5  10  15 Inf];

%%
% A time gated spectral measurement of the Non-HT Data field is used for
% the transmitter spectrum emission mask test [ <#11 3> ]. The Non-HT Data
% field of each packet is extracted from the upsampled |txWaveform| using
% the start index of each packet. The extracted Non-HT Data fields are
% concatenated in preparation for measurement.

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgNHT);
startIdx = osf*(ind.NonHTData(1)-1)+1;   % Upsampled start of Non-HT Data
endIdx = osf*ind.NonHTData(2);           % Upsampled end of Non-HT Data
idleNSamps = osf*idleTime/(1/fs);        % Upsampled idle time samples
perPktLength = endIdx + idleNSamps;

idx = zeros(endIdx-startIdx+1, numPackets);
for i = 1:numPackets
    % Start of packet in txWaveform, accounting for the filter delay
    pktOffset = (i-1)*perPktLength+filterLen/2;
    % Indices of NonHT Data in txWaveform
    idx(:,i) = pktOffset+(startIdx:endIdx);
end
% Select the Data field for the individual packets
gatedNHTDataTx = txWaveform(idx(:),:);

%%
% The plot generated by the helper function |helperSpectralMaskTest|
% overlays the required spectral mask with the measured power spectral
% density. It checks the transmitted PSD levels to be within the
% specified mask levels and displays a pass/fail status after the test.

% Evaluate the PSD and check for compliance
helperSpectralMaskTest(gatedNHTDataTx, fs, osf, dBrLimits, fLimits);

% Restore default stream
rng(s);

%% Conclusion and Further Exploration
% The transmit spectral mask for Class A Stations at the 5.85-5.925 GHz
% bands for a 10 MHz channel spacing is shown in this example. It is also
% shown how the peak spectral density of the transmitted signal falls
% within the spectral mask to satisfy regulatory restrictions. A similar
% result can be generated for the 5 MHz channel spacing.
%
% The high power amplifier model and the spectral filtering affect the
% out-of-band emissions in the spectral mask plot. For different station
% classes with higher relative dB values, try using different filters or
% filter lengths and/or increase the backoff for lower emissions.
%
% For information on other transmitter measurements like modulation
% accuracy and spectral flatness, refer to the following example:
%
% * <VHTTransmitterMeasurementsExample.html 802.11ac(TM) Transmitter
% Modulation Accuracy and Spectral Emission Testing>

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperSpectralMaskTest.m') helperSpectralMaskTest.m>
% * <matlab:edit('helperSampleRate.m') helperSampleRate.m>

%% Selected Bibliography
% # IEEE Std 802.11-2012: IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements, Part 11: Wireless LAN
% Medium Access Control (MAC) and Physical Layer (PHY) Specifications,
% IEEE, New York, NY, USA, 1999-2013.
% # IEEE Std 802.11p-2010: IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements, Part 11: Wireless LAN
% Medium Access Control (MAC) and Physical Layer (PHY) Specifications,
% Amendment 6: Wireless Access in Vehicular Environments, IEEE, New York,
% NY, USA, 2010.
% # Archambault, Jerry, and Shravan Surineni. "IEEE 802.11 spectral
% measurements using vector signal analyzers." RF Design 27.6 (2004):
% 38-49.

displayEndOfDemoMessage(mfilename)