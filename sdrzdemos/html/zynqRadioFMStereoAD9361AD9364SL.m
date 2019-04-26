%% FM Stereo Receiver Using Analog Devices AD9361/AD9364
% This model shows how to use the Xilinx(R) Zynq-Based Radio Support
% Package with Simulink(R) to build an FM Stereo receiver.  The receiver
% receives and demodulates the FM signal transmitted by the stereo FM
% broadcast radio.
%
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting
% Started> documentation for details on configuring your host computer to
% work with the Support Package for Xilinx(R) Zynq-Based Radio.

% Copyright 2014-2015 The MathWorks, Inc.

%% Introduction
% Xilinx Zynq-Based SDR hardware along with Analog Devices(TM)
% AD9361/AD9364 is required for running the model. The model configures the
% SDR hardware to receive stereo FM signal over the air. It plays the left
% and right channels.

%% Structure of the Example
%
% The top level structure of the FM Mono receiver model is shown in the
% figure below:
modelName = mfilename;
open_system(modelName);
set_param(modelName, 'SimulationCommand', 'update');

%%
% In this example, the radio baseband data is received from the SDR
% Hardware via a Zynq(R) SDR Receiver block
% <matlab:sdrzdoc('comminternalsdrrxzc706fmc23sl') ZC706 and FMCOMMS2/3/4
% Receiver>. The Zynq SDR receiver block is configured at a sampling rate
% of 960 kHz using the _Baseband sample rate_ parameter in the Zynq SDR
% Receiver block. 

% The |FM Broadcast Demodulator Baseband| block converts the sampling 
% rate of 960 kHz to 48 kHz, a native sampling rate for your host 
% computer's audio device. According to the FM broadcast standard in the 
% United States,  the deemphasis lowpass filter time constant is set 
% to 75 microseconds.

close_system(modelName,0);
%%
% To perform stereo decoding, the |FM Broadcast Demodulator Baseband| block
% uses a peaking filter which picks out the 19 kHz pilot tone from which
% the 38 kHz carrier is created. Using the obtained carrier signal, the
% block downconverts the L-R signal, centered at 38 kHz, to baseband.
% Afterwards, the L-R and L+R signals pass through a 75 microsecond
% deemphasis filter . The |FM Broadcast Demodulator Baseband| block
% separates the L and R signals and converts them to the 48 kHz audio
% signal
%
%% Running the Example
%
% The SDR hardware is configured to use input ports for the RF center
% frequency. Using input ports for these parameters allows them to be
% updated while the model is running. The center frequency can be
% configured to a local FM radio station using the *Slider Gain* block. Set
% the center frequency to a local FM radio station and start the model.
% This tunes the SDR device to the local FM station and you can listen to
% the radio station through the audio device connected to your computer.
% Change the center frequency to listen to a different station.
%
% If you hear some dropouts or delay in the sound, run the model in
% Accelerator mode. From the model menu, select Simulation->Accelerator,
% then click the run button.
%

%% Exploring the Example
% If you have your own FM transmitter that can transmit .wma files, you can
% duplicate the test that shows the channel separation result above. Load
% the |sdrzFMStereoTestSignal.wma| file into your transmitter. The
% channel separation can be easily observed from the Spectrum Scope block and
% heard from the audio device. 
%
%% Troubleshooting the Example
%
% If the received signal is very weak, you can try increasing the receiver
% gain by changing the _Source of gain_ variable to _Input port_ for the
% manual gain control mode or by changing it to 'AGC Fast Attack' or 'AGC
% Slow Attack'.
%
% If you run the example as described but fail to see a signal like the one
% shown (e.g. you only receive noise or the spectrum display is never
% updated), see <matlab:sdrzdoc('sdrz_troubleshoot') Xilinx FPGA-Based
% Radio Processing Errors and Fixes>.
%
%% References
% * <http://en.wikipedia.org/wiki/FM_broadcasting FM broadcasting on Wikipedia(R)>
% 
