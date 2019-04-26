%% Walkie-Talkie Transmitter with Analog Devices FMCOMMS1
%
% This example shows how to use the Xilinx(R) Zynq-Based Radio Support
% Package with Simulink(R) to implement a walkie-talkie transmitter. The
% transmitted signal can be received by a compatible commercial
% walkie-talkie. Alternatively, the signal can be received by the companion
% <matlab:showdemo('sdrzWalkieTalkieRx') Walkie-Talkie Receiver with
% Analog Devices(TM) FMCOMMS1> example if you have a second SDR platform.
%
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting
% Started> documentation for details on configuring your host computer to
% work with the Support Package for Xilinx(R) Zynq-Based Radio.

% Copyright 2014 The MathWorks, Inc.


%% Introduction
%
% Walkie-talkies provide a subscription-free method of communicating over
% short distances. There are a number of different standards used around
% the world. This example uses Simulink to implement two of these
% standards: *Family Radio Service* and *Personal Mobile Radio 446*.
%
% * *Family Radio Service (FRS):* Operates on 14 channels at frequencies
% around 462 MHz and 467 MHz. The channel spacing is 25 kHz. FRS radios are
% commonly found in the North and South America. More details on FRS can be
% found in [ <#21 1> ].
% * *Personal Mobile Radio 446 (PMR446):* Operates on 8 channels around 446
% MHz. The channel spacing is 12.5 kHz. PMR446 radios are commonly found in
% Europe. More details on PMR446 can be found in [ <#21 2> ].
%
% Both FRS and PMR446 use analog frequency modulation (FM) with a maximum
% frequency deviation of +-2.5 kHz. They also use a *Continuous Tone-Code
% Squelch System (CTCSS)* to filter unwanted signals at the receiver. CTCSS
% is implemented in this example.
%
% This example allows the transmitted audio to be a continuous tone, a
% chirp signal or voice from an audio file. The audio sample rate is always
% 8 kHz.


%% Setup
%
% Before running the example, ensure you have performed the following
% steps:
%
% 1. Configure your host computer to work with the Support Package for
% Xilinx(R) Zynq-Based Radio. See <matlab:sdrzdoc('sdrzspsetup')
% Getting Started> for help.
%
% * Some additional steps may be required if you want to run two radios
% from a single host computer. See
% <matlab:sdrzdoc('sdrz_tworadios') Setup for Two Radios - One
% Host> for help.
%
% 2. Ensure that you have a suitable receiver. This example is designed to
% work in conjunction with any of the following possible receivers:
%
% * The <matlab:showdemo('sdrzWalkieTalkieReceiver') Walkie-Talkie
% Receiver with Analog Devices FMCOMMS1> MATLAB(R) example
% * The <matlab:showdemo('sdrzWalkieTalkieRx') Walkie-Talkie Receiver with
% Analog Devices FMCOMMS1> Simulink example
% * A commercial FRS/PMR446 radio
%
% 3. Ensure that your receiver is set to the same protocol, channel and
% CTCSS code.


%% Running the Example
%
% Before running the example, make sure that walkie-talkie protocol is set
% to a standard that is legal for use in your location.
%
% Once both the transmitter and receiver are active, you should hear the
% transmitted audio on your receiver. By default, a voice signal is
% transmitted.
%
% You can configure various settings that control the model using the
% *Model Parameters* block. Except for the walkie-talkie protocol, all of
% the parameters can be changed while the model is running.


%% Transmitter Design: System Architecture
%
% The top level of the model is shown below.
modelName = mfilename;
open_system(modelName);
set_param(modelName, 'SimulationCommand', 'update')
%%
% The system has been split up into a number of subsystems, which each
% perform a specific task:
%
% # *Audio Generator:* Generate the audio signal to be transmitted
% # *CTCSS Generator:* Generate the CTCSS tone
% # *Software Interpolator:* Upsample the 8 kHz audio signal to 64 kHz
% # *FM Baseband Modulator:* FM modulate the audio signal
% # *SDR Transmitter:* Send the modulated signal to the SDR hardware for
% further upsampling and transmission
%
% Each subsystem is explored in more detail below.
%


%%
% *Audio Generator*
%
open_system([modelName '/AudioSources']);
%%
% The *Audio Generator* subsystem is where the audio signal that is audible
% at the receiver is generated. A simple Multiport Switch block allows the
% audio source to be switched during model execution. The audio signal can
% be a single tone, a chirp or audio played from a multimedia file.
%


%%
% *CTCSS Generator*
%
close_system([modelName '/AudioSources']);
open_system([modelName '/CTCSSSource']);
%%
% The *CTCSS Generator* subsystem is where the CTCSS tone is generated. The
% Continuous Tone-Coded Squelch System (CTCSS)[ <#21 3> ] allows receivers
% to filter out received signals that did not originate from the target
% transmitter. The transmitter generates a tone between 67 Hz and 250 Hz
% and transmits it along with the audio signal. The receiver contains logic
% to detect this tone, and acknowledges a message if the detected tone
% matches the code setting on the receiver. The receiver filters out the
% tone so that the user does not hear it.
%
% The CTCSS tone generator generates a continuous phase sine wave with a
% frequency corresponding to the selected private code. The amplitude of
% the tone is usually 10%-15% of the maximum amplitude of the modulating
% signal. Note that because the maximum amplitude of all the source signals
% is 1, the default amplitude of 0.15 for the CTCSS tone corresponds to 15%
% of the modulating signal's maximum amplitude.

%%
% *Software Interpolator*
%
close_system([modelName '/CTCSSSource']);
open_system([modelName '/SoftwareInterpolator']);
%%
% The *Software Interpolator* subsystem is where the 8 kHz audio signal is
% upsampled to the radio baseband rate of 64 kHz. It uses an FIR
% Interpolator block, the filter coefficients of which are generated in the
% helper function <matlab:edit('sdrzWalkieTalkieTxHelper_SimParams.m')
% sdrzWalkieTalkieTxHelper_SimParams.m>. 

%%
% *FM Baseband Modulator*
%
close_system([modelName '/SoftwareInterpolator']);
open_system([modelName '/FMBasebandModulator']);
%%
% The *FM Baseband Modulator* subsystem is where the upsampled audio signal
% is frequency modulated. This example implements the FM modulator using a
% simple digital IIR filter as an integrator. The frequency sensitivity
% gain is used to control the modulation. It is related to the frequency
% deviation by the formula:
%
% _frequencySensitivityGain = frequencyDeviation * (2*pi*Ts) / A_
%
% where _peakFrequencyDeviation_ is 2.5 kHz, _Ts_ is the sampling period of
% the SDR transmitter, and _A_ represents the maximum amplitude of the
% modulating signal i.e. the audio signal. This example assumes the
% generated audio is normalized and therefore has a peak amplitude of 1.
%
% See [ <#21 4> ] for more information on frequency modulation.

%%
% *SDR Transmitter*
%
close_system([modelName '/FMBasebandModulator']);
open_system([modelName '/SDRTransmitter']);
%%
% The *SDR Transmitter* subsystem is where the radio baseband data is
% passed to the SDR hardware via an
% <matlab:sdrzdoc('comminternalsdrtxzc706fmc1txsl') ZC706 and Analog Devices FMCOMMS1
% Transmitter> block. The SDR hardware is configured to use input ports for
% the RF center frequency and the intermediate frequency. Using an input
% port for these parameters allows them to be updated while the model is
% running. The desired center frequency is dependent on the channel
% selected in the *Model Parameters* block. The intermediate frequency is
% defined in the helper function
% <matlab:edit('sdrzWalkieTalkieTxHelper_SimParams.m')
% sdrzWalkieTalkieTxHelper_SimParams.m> to be 10 MHz.


%% Transmitter Design: System Sample Rates
%
% The system has three different sample rates:
%
% # The audio sample rate, *8 kHz*
% # The SDR hardware baseband rate, *64 kHz*
% # The SDR hardware DAC rate, *32.768 MHz*
%
% The upsample by 8 from 8 kHz to 64 kHz is done in software. The upsample
% to 64 kHz is necessary for two reasons:
%
% # By Carson's rule, the approximate passband bandwidth of the desired FM
% signal is 2*(frequencyDeviation + peakAudioFrequency). In this example,
% that equates to a signal bandwidth of _2*(2.5e3 + 4e3) = 13 kHz_. This
% means we need to use a sample rate greater than 13 kHz.
% # The sample rates in software have no relation to the real world
% transmission rates. The actually transmission rate is determined entirely
% by the SDR hardware DAC rate and the SDR hardware interpolation factor.
% To make the software sample rates meaningful, the sample rates at the
% software/hardware interface must match. For the FMCOMMS1, the lowest
% possible baseband sample rate that is a multiple of 8 kHz that is greater
% than 13 kHz is 64 kHz. Using an integer upsample factor means a smaller
% interpolation filter (fewer taps) can be used.
%

%% Things to Try
%
% To modify the example, double click on the *Model Parameters* block. Most
% of the parameters can be altered while the model is running. Some
% possible modifications include:
%
% * Try changing the channel.
% * Try changing the CTCSS code. Note that the receiver will not play the
% transmission out loud unless it has the same CTCSS code, or it has CTCSS
% disabled by setting the code to 0.
% * Try changing the audio source. You can send a voice recording, a single
% tone or a chirp signal.
% * Try replacing the _From Multimedia File_ block in the *Audio Generator*
% subsystem with a _From Audio Device_ block from the DSP System Toolbox(TM).
% This would allow you to transmit audio captured from a
% microphone in real time.


%% Alternative Implementations
%
% This example implements a walkie-talkie transmitter in Simulink. You can
% also view the equivalent system implemented using MATLAB in the
% <matlab:showdemo('sdrzWalkieTalkieTransmitter') Walkie-Talkie Transmitter
% with Analog Devices FMCOMMS1 using MATLAB> example.


%% Troubleshooting the Example
%
% If you cannot hear the transmitted signal on your receiver, try the
% following:
%
% * Make sure that the transmitter and the receiver are set to the same
% channel and CTCSS code.
% * Disable CTCSS on the receiver by setting the code to 0. Note that
% codes higher than 38 use Digital Coded Squelch, which is not implemented
% in this example.
%
% General tips for troubleshooting SDR hardware can be found in
% <matlab:sdrzdoc('sdrz_troubleshoot') Xilinx Zynq-Based Radio Processing
% Errors and Fixes>.


%% List of Example Helper Files
%
% This example uses the following helper files
%
% * <matlab:edit('sdrzWalkieTalkieTxHelper_SimParams.m')
% sdrzWalkieTalkieTxHelper_SimParams.m>: returns a structure of variables
% used to control constant parameters in the model.
% * <matlab:edit('sdrzWalkieTalkieTxHelper_ModelParamsMask.m')
% sdrzWalkieTalkieTxHelper_ModelParamsMask.m>: controls the *Model
% Parameters* mask and updates the model depending on the parameters
% supplied in the mask.
% * <matlab:edit('sdrzWalkieTalkieHelper_CTCSSCode2Tone.m')
% sdrzWalkieTalkieHelper_CTCSSCode2Tone.m>: converts a CTCSS code to a
% frequency.
% * <matlab:edit('sdrzWalkieTalkieHelper_Channel2Frequency.m')
% sdrzWalkieTalkieHelper_Channel2Frequency.m>: converts a channel number to
% an RF frequency. 
% * sdrzWalkieTalkieHelper_voice.wav: the audio file used when the
% audio source is set to _'Audio file'_.

%% References
%
% # <http://en.wikipedia.org/wiki/Family_Radio_Service Family Radio
% Service> on Wikipedia(R)
% # <http://en.wikipedia.org/wiki/PMR446 PMR446> on Wikipedia
% # <http://en.wikipedia.org/wiki/Continuous_Tone-Coded_Squelch_System
% Continuous Tone-Coded Squelch System> on Wikipedia
% # <http://en.wikipedia.org/wiki/Frequency_modulation Frequency Modulation>
% on Wikipedia

%%
close_system([modelName '/SDRTransmitter']);
close_system(modelName, 0);