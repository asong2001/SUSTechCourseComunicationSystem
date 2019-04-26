%% FM Monophonic Receiver Using Analog Devices AD9361/AD9364
% This model shows how to use the Xilinx(R) Zynq-Based Radio Support
% Package with Simulink(R) to build an FM mono receiver. The receiver
% receives and demodulates the FM signal transmitted by the FM broadcast radio.
%
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting
% Started> documentation for details on configuring your host computer to
% work with the Support Package for Xilinx(R) Zynq-Based Radio.

% Copyright 2014-2015 The MathWorks, Inc.

%% Introduction
%
% Xilinx Zynq-Based SDR hardware along with Analog Devices(TM)
% AD9361/AD9364 is required for running the model. The model configures
% the SDR hardware to receive an FM signal over the air. The FM Mono
% receiver performs FM demodulation, deemphasis filter and rate conversion.
%
%% Setup
%
% Before running the example, ensure you have performed the following
% steps:
%
% Configure your host computer to work with the Support Package for
% Xilinx(R) Zynq-Based Radio. See <matlab:sdrzdoc('sdrzspsetup')
% Getting Started> for help.
%
% * Some additional steps may be required if you want to run two radios
% from a single host computer. See
% <matlab:sdrzdoc('sdrz_tworadios') Setup for Two Radios - One
% Host> for help.
%
%% Structure of the Example
%
% The top level structure of the FM Mono receiver model is shown in the figure below:
modelName = mfilename;
open_system(modelName);
set_param(modelName, 'SimulationCommand', 'update');

%%
% In this example, the radio baseband data is received from the SDR Hardware via
% a Zynq(R) SDR Receiver block  <matlab:sdrzdoc('comminternalsdrrxzc706fmc23sl') ZC706
% and FMCOMMS2/3/4 Receiver>. The Zynq SDR receiver block is configured at a sampling rate of 
% 960 kHz using the _Baseband sample rate_ parameter in the Zynq SDR Receiver block. The
% deemphasis filter in the *FM Receiver* subsystem is set to 75 microseconds [ <#11 1> ]
% and compensates for the preemphasis filter at the transmitter.  The
% frequency response table is given below.
%
% <<zynqRadioFMMonoAD9361AD9364_deemp.png>>
%
% A resampler converts the sampling rate from 960 kHz to 48 kHz, a native
% sampling rate for the audio device. The resampling filter removes the 19
% kHz stereo pilot tone.
%

%% Running the Example
%
% The SDR hardware is configured to use input ports for the RF center frequency.
% Using input ports for these parameters allows them to be
% updated while the model is running. The center frequency can be
% configured to a local FM radio station using the *Slider Gain* block.
%
% Set the center frequency to a local FM radio station and start the model. 
% This tunes the SDR device to the local FM station and you can listen to the 
% radio station through the audio device connected to your computer.  Change the
% center frequency to listen to a different station.
%
%%
close_system(modelName,0);
%% Alternative Implementations
%
% The example describes the Simulink model of a receiver for
% demodulating and processing an FM Mono signal transmitted by the FM broadcast
% radio. The example uses Zynq-based Hardware support package to build a monophonic 
% FM receiver using SDR Hardware using Analog Devices AD9361/AD9364.
%
% You can also view a MATLAB(R) implementation of the example in
% <matlab:showdemo('zynqRadioFMMonoReceiverAD9361AD9364ML') FM Monophonic Receiver Using Analog Devices AD9361/AD9364 using MATLAB> 
%
%% Troubleshooting the Example
%
% If the received signal is very weak, you can try increasing the receiver
% gain by changing the _Source of gain_ variable to _Input port_ for the manual gain control mode
% or by changing it to 'AGC Fast Attack' or 'AGC Slow Attack'.
%
% General tips for troubleshooting SDR hardware can be found in
% <matlab:sdrzdoc('sdrz_troubleshoot') Xilinx Zynq-Based Radio Processing
% Errors and Fixes>.
%
%% References
% * <http://en.wikipedia.org/wiki/FM_broadcasting FM broadcasting on Wikipedia(R)>
%
