%% QPSK Transmitter Using Analog Devices AD9361/AD9364
%
% This example shows how to use the Xilinx(R) Zynq-Based Radio Support
% Package with Simulink(R) to implement a QPSK transmitter. The SDR device
% in this model will continuously transmit indexed 'Hello world' messages
% that are QPSK modulated onto a carrier with a specified center frequency.
% You can demodulate the transmitted message using the
% <matlab:showdemo('zynqRadioQPSKRxAD9361AD9364SL') QPSK Receiver Using Analog Devices(TM) AD9361/AD9364> model if
% you have a second SDR platform.
%
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting Started>
% documentation for details on configuring your host computer to 
% work with the Support Package for Xilinx(R) Zynq-Based Radio.

% Copyright 2014 The MathWorks, Inc.

%% Introduction
% 
% This example transmits a QPSK signal over the air using SDR hardware. It
% has two main objectives:
%
% * Implement a prototype QPSK-based transmitter in Simulink(R) using
% Simulink blocks from the  Xilinx(R) Zynq-Based Radio Support Package.
%
% * Illustrate the use of key Communications System Toolbox(TM) Simulink
% blocks for QPSK system design.


%% Setup
% 
% Before running the example, ensure you have performed the following
% steps:
%
% 1. Configure your host computer to work with the Support Package for
% Xilinx(R) Zynq-Based Radio. See <matlab:sdrzdoc('sdrzspsetup') Getting Started> for help.
%
% * Some additional steps may be required if you want to run two radios
% from a single host computer. See
% <matlab:sdrzdoc('sdrz_tworadios') Setup for Two Radios-One Host>
% for help.
%
% 2. Ensure that you have a suitable receiver. This example
% is designed to work in conjunction with any of the following possible
% receiver examples:
%
% * The  <matlab:showdemo('zynqRadioQPSKReceiverAD9361AD9364ML') QPSK Receiver Using Analog Devices AD9361/AD9364>
% MATLAB(R) example
% * The <matlab:showdemo('zynqRadioQPSKRxAD9361AD9364SL') QPSK Receiver Using Analog Devices AD9361/AD9364> Simulink example
% * The <matlab:showdemo('zynqRadioQPSKRxFPGAAD9361AD9364SL') Targeting HDL-Optimized QPSK Receiver Using Analog Devices AD9361/AD9364> Simulink example

%% Running the Example
%
% Start the transmitter and then your companion receiver. Once both are
% running, you should see "Hello world" messages in the MATLAB command
% window where the receiver is running.

%% Transmitter Design: System Architecture
%
% The top-level structure of the model is shown below.
modelname = 'zynqRadioQPSKTxAD9361AD9364SL';
open_system(modelname);
set_param(modelname, 'SimulationCommand', 'update')

%%
% The system performs performs four major processes: 
%
% # Bit generation
% # Baseband modulation
% # Pulse shaping and upsampling
% # Sending baseband data to SDR hardware
% 
% Each process is explored in more detail in the following sections. The
% result of this processing is a single channel of QPSK data for 
% transmission by the Zynq(R) SDR Transmitter block.
%
% *Bit Generation*
%
% The *Bit Generation* subsystem uses a MATLAB workspace variable as the
% payload of a frame, as shown in the figure below. 
%
open_system([modelname '/Bit Generation']);


%%
% Each frame contains 200 bits. The first 26 bits are a frame header, and
% the remaining 174 bits represent a data payload.
%
% * The 26 header bits result in a
% 13-symbol Barker code to use as a preamble. The preamble is used to
%  aid in overcoming channel impairments in the receiver. 
% * The first 105 bits of the payload
% correspond to the ASCII representation of 'Hello world ###', where '###'
% is a repeating sequence of '001', '002', '003',..., '099'.  
% * The remaining
% payload bits are random. 
%
% The payload is scrambled to guarantee a
% balanced distribution of zeros and ones for the timing recovery operation
% in the receiver.
%
% *Baseband Modulation*
%
% The *QPSK Modulator Baseband* block modulates pairs of bits from the
% output of the *Bit Generation* subsystem to QPSK constellation points
% using Gray mapping. Each QPSK symbol is represented by one complex
% sample.
%
% *Pulse Shaping and Upsampling*
%
% The *Raised Cosine Transmit Filter* block performs root raised cosine
% pulse shaping with a roll off factor of 0.5. It also upsamples the
% baseband signal by a factor of 4.
%
% *Sending Baseband Data to SDR Hardware*
%
% The <matlab:sdrzdoc('comminternalsdrtxzc706fmc23sl') ZC706 and Analog Devices
% FMCOMMS2/3/4 Transmitter> block sends baseband data to the SDR hardware over ethernet.
% The FPGA sends the baseband data to match the AD9361/AD9364 baseband sample rate, at
% which point the AD9361/AD9364 further upsamples the signal to RF and transmits
% it over the air. It is important to note that the real world rate at
% which the model runs is determined by the *Baseband sample rate* in 
% *Zynq SDR Transmitter* block, and *not* by the simulation sample time. 
close_system(modelname, 0);

%% Alternative Implementations
%
% This example describes the Simulink implementation of a QPSK transmitter
% with Zynq SDR platform and Analog Devices AD9361/AD9364. You can also view a MATLAB implementation of this
% example in <matlab:showdemo('zynqRadioQPSKTransmitterAD9361AD9364ML') QPSK Transmitter with
% Analog Devices AD9361/AD9364 using MATLAB>.
%
% You can also explore a non-hardware QPSK transmitter and receiver example
% that models a general wireless communication system using an AWGN channel
% and simulated channel impairments with the
% <matlab:showdemo('commqpsktxrx') QPSK Transmitter and Receiver> example.


%% Troubleshooting the Example
% 
% If you run the example and you get the message |WARNING: SDR hardware Tx
% data buffer underflow!| in the command window, then the simulation ran
% slower than real time. You can try using
% <matlab:sdrzdoc('sdrz_burstmode') burst mode>.
% 
% If you still fail to get the example to work, see
% <matlab:sdrzdoc('sdrz_troubleshoot') Xilinx FPGA-Based Radio Processing
% Errors and Fixes>.


%% List of Example Helper Files
%
% This example uses the following helper files:
%
% * <matlab:edit('zynqRadioQPSKTxAD9361AD9364SL_init.m') zynqRadioQPSKTxAD9361AD9364SL_init.m>: returns a
% structure of variables used to control constant parameters the model.
