%% Frequency Offset Calibration with Analog Devices FMCOMMS1
%
% This example shows how to use the Xilinx(R) Zynq-Based Radio Support
% Package with Simulink(R) to determine the frequency offset between SDR
% devices. The transmitter sends a 10 kHz sine wave with the
% <matlab:sdrzfreqcalib Frequency Offset Calibration (Tx) with SDR
% Hardware> model. The receiver receives the signal, calculates the
% frequency offset and displays the offset in the <matlab:sdrzfreqcalib_rx
% Frequency Offset Calibration (Rx) with Analog Devices(TM) FMCOMMS1> model.
%
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting
% Started> documentation for details on configuring your host computer to
% work with the Support Package for Xilinx(R) Zynq-Based Radio.

% Copyright 2014 The MathWorks, Inc.


%% Introduction
%
% This example uses a matched pair of models to determine the frequency
% offset between two SDR devices:
%
% * The transmit model is <matlab:sdrzfreqcalib sdrzfreqcalib>
% * The receive model is <matlab:sdrzfreqcalib_rx sdrzfreqcalib_rx>
%
% The transmitter sends a 10 kHz tone. The receiver detects the transmitted
% tone using an FFT-based detection method. The offset between the
% transmitted 10 kHz tone and the received tone can then be calculated and
% used to compensate for the offset at the receiver. The pair of models
% provides the following information:
%
% * A quantitative value of the frequency offset
% * A graphical view of the spur-free dynamic range of the receiver
% * A graphical view of the qualitative SNR level of the received signal


%% Setup
% 
% Before running the example, make sure you have performed the following
% steps:
%
% 1. Configure your host computer to work with the Support Package for
% Xilinx(R) Zynq-Based Radio. See <matlab:sdrzdoc('sdrzspsetup')
% Getting Started> for help.
%
% * Some additional steps may be required if you want to run two radios
% from a single host computer. See
% <matlab:sdrzdoc('sdrz_tworadios') Setup for Two Radios-One Host>
% for help.
%
% 2. Make sure that you have both the transmitter model
% <matlab:sdrzfreqcalib sdrzfreqcalib> and the receiver model
% <matlab:sdrzfreqcalib_rx sdrzfreqcalib_rx> open, with each configured to
% run on its own SDR hardware in its own instance of Simulink.


%% Running the Example
%
% Start the transmitter model running, and then start the receiver model.
%
% The calculated frequency offset is shown by the *Frequency Offset
% Display* block in the receiver model. The *Spectrum Analyzer* block in
% the *Receiver* subsystem shows the spectrum of the received signal. A
% sample spectrum is shown below.
%
% <<sdrzfreqcalibspectrum.png>>
%
% In this case, the frequency with the maximum received signal power is at
% about 9.37 kHz. Since the transmitter is sending a tone at 10 kHz, this
% means the frequency offset is about 630 Hz. The spurious free dynamic
% range of the signal is about 49 dB.
%
% To compensate for a transmitter/receiver frequency offset, add the
% displayed frequency offset value to the _Intermediate frequency_
% parameter of the <matlab:sdrzdoc('comminternalsdrtxzc706fmc1txsl') ZC706 and Analog Devices FMCOMMS1
% Transmitter> block. Be sure to use
% the sign of the offset in your addition. Rerun the receiver with the
% applied frequency offset compensation. The calculated offset frequency
% displayed should now be close to zero, and the peak in the spectrum
% should be close to 10 kHz.
%
% The offset correction is applied to the _Intermediate frequency_
% parameter instead of the _Center frequency_ because the _Intermediate
% frequency_ provides a much finer tuning resolution. 
%
% It is important to note that the frequency offset value is only valid for
% the center frequency used to run the calibration.


%% Transmitter Design: System Architecture
%
% The following figure shows the transmitter model:
modelname1 = 'sdrzfreqcalib';
open_system(modelname1);
set_param(modelname1, 'SimulationCommand', 'update')
%%
% The transmitter sends a 10 kHz tone at a default center frequency of
% 2.415 GHz i.e. the tone is transmitted at 2.415 GHz + 10 kHz. 

%% Receiver Design: System Architecture
%
% The following figure shows the receiver model:
modelname2 = 'sdrzfreqcalib_rx';
close_system(modelname1, 0);
open_system(modelname2);
close_system([modelname2 '/Receiver/Spectrum Analyzer']); % Don't publish the empty scope
set_param(modelname2, 'SimulationCommand', 'update')
%%
% The following figure shows the detailed structure of the *Receiver*
% subsystem:
open_system([modelname2 '/Receiver']);
%%
% * The *Spectrum Analyzer* block computes and displays the power
% spectral density of the received signal.
% * The *Find Peak Frequency* subsystem uses an FFT to find the frequency
% with the maximum power in the received signal.
%
% The *Find Peak Frequency* subsystem finds the frequency with the maximum
% power in the received signal. The subsystem is shown below:
close_system([modelname2 '/Receiver']);
open_system([modelname2 '/Receiver/Find Peak Frequency'],'force');
%%
% In this subsystem, the *Periodogram* block returns the PSD estimate of
% the received signal. The *Probe* block finds the frame size and the frame
% sample time. With this information, the MATLAB(R) function block
% *findpeakfreq* finds the index of the maximum amplitude across the
% frequency band and converts the index to a frequency value according to
% the following calculation:
%
%  Foffset = IndexofMaxAmplitude * FrameSize / (FFTLength * FrameSampleTime)
%%
% Note that the FFT performed by the *Periodogram* uses 4096 samples. This
% means that the frequency offset calculated is limited to a resolution of
% 48 Hz.


%% Alternative Implementations
%
% This example describes the Simulink implementation of a pair of models
% for performing frequency offset calibration between two SDR devices. You
% can also view a MATLAB implementation of these models in
% <matlab:showdemo('sdrzFrequencyCalibrationTransmitter') Frequency Offset
% Calibration Transmitter with Analog Devices FMCOMMS1 using MATLAB> and
% <matlab:showdemo('sdrzFrequencyCalibrationReceiver') Frequency Offset
% Calibration Receiver with Analog Devices FMCOMMS1 using MATLAB>.

%% Troubleshooting the Example
%
% If the received signal is very weak, you can try increasing the receiver
% gain by increasing the _Gain_ value in the *Analog Devices FMCOMMS1
% Receiver* block.
%
% If you run the example as described but fail to see a signal like the one
% shown (e.g. you only receive noise or the spectrum display is never
% updated), see <matlab:sdrzdoc('sdrz_troubleshoot') Xilinx FPGA-Based
% Radio Processing Errors and Fixes>.

close_system([modelname2 '/Receiver/Find Peak Frequency']);
close_system(modelname2, 0);
