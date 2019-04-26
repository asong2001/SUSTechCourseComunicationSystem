%% Targeting HDL Optimized QPSK Transmitter with Analog Devices FMCOMMS1
%
% This example shows how to model an HDL-optimized QPSK transmitter and
% prototype it on the SDR hardware using the HDL Coder(TM) workflow advisor.
% The SDR device in this model will keep transmitting indexed 'Hello world'
% messages at its specified center frequency. For a full description of
% targeting, see <matlab:sdrzdoc('sdrz_targeting') Targeting
% Overview> in the documentation.
%
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting Started>
% documentation for details on configuring your host computer to work with
% the Support Package for Xilinx(R) Zynq-Based Radio. Additionally,
% targeting requires HDL Coder and Xilinx Vivado(R) 2015.2.1, which must be
% on the system path. You can add Vivado to the path when you call
% <matlab:sdrzdoc('sdrz_settoolpath') setuptool>. See the
% <matlab:sdrzdoc('sdrz_targeting') targeting requirements> for more
% details.

% Copyright 2014-2015 The MathWorks, Inc.


%% Introduction
%
% The <matlab:showdemo('commqpsktxhdl') commqpsktxhdl> example demonstrates
% building a practical HDL-optimized QPSK transmitter, which includes QPSK
% modulation and root raised cosine pulse shaping with an oversampling
% ratio of 4. This example is a combination of the
% <matlab:showdemo('sdrzqpsktx') sdrzqpsktx> example and the commqpsktxhdl
% example. It takes the next steps to generate HDL code from the
% HDL-optimized transmitter and prototype it on the FPGA of the SDR
% hardware. Generation of bit frames is done in Simulink(R) before being
% packed for transmission over Ethernet to the SDR hardware. The raw data
% bits are modulated and the resulting symbols are pulse shaped in the
% programmable logic before being upsampled and sent to the RF front end.
%
% Targeting typically involves a pair of models:
%
% # The targeting model, used to generate the HDL code and integrate the
% DUT (Design Under Test) into the transmit path implemented on the FPGA of
% the SDR hardware
% # The retargeted model, which uses the SDR hardware with the included DUT
% performing a portion of the baseband processing
%
% Note that there are two methods of running this example, based on the
% targeting/retargeted pair of models:
%
% # Run the <matlab:sdrzQPSKTxFPGA sdrzQPSKTxFPGA> targeting model
% directly. See the <#3 Running the Example: Simulated> section.
% # Generate a targeted bitstream and use the retargeted model
% <matlab:sdrzQPSKTxFPGAretarget sdrzQPSKTxFPGAretarget>. See the <#4
% Running the Example: Targeted> section.


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
% * The <matlab:showdemo('sdrzqpskrx') QPSK Receiver with Analog Devices(TM) FMCOMMS1>
% Simulink example
% * The <matlab:showdemo('sdrzQPSKReceiver') QPSK Receiver with Analog Devices FMCOMMS1> MATLAB(R) example
% * The <matlab:showdemo('sdrzQPSKRxFPGA') Targeting HDL-Optimized QPSK
% Receiver with Analog Devices FMCOMMS1> Simulink example


%% Running the Example: Simulated
%
% To run the example directly, open <matlab:sdrzQPSKTxFPGA sdrzQPSKTxFPGA>
% and start the simulation. Note that this simulates the entire baseband
% design in a similar way to the <matlab:showdemo('sdrzqpsktx') QPSK
% Transmitter with Analog Devices FMCOMMS1> Simulink example i.e. no HDL code is
% generated, and all of the baseband processing is done in Simulink.
%
% The main benefit of running the example in this way is that you can
% verify the conversion of your design to fixed point HDL works with real
% world data. The main disadvantage is that the HDL optimized sections of
% the model will run *significantly* slower than their floating point,
% frame based equivalents. This example uses
% <matlab:sdrzdoc('sdrz_burstmode') burst mode> to ensure large enough
% blocks of data are received consecutively to verify the receiver
% functionality.


%% Running the Example: Targeted
%
% The <matlab:sdrzQPSKTxFPGA sdrzQPSKTxFPGA> model can be used to generate
% a targeted SD card image for the SDR hardware that includes the *HDLTx*
% subsystem as a DUT. The DUT is automatically placed in the SDR transmit
% chain that is implemented on the FPGA.
%
% To generate the new SD card image, perform the following steps:
%
% # Make sure you have added Vivado to the system path. See
% <matlab:sdrzdoc('sdrz_settoolpath') setuptool>.
% # Open the HDL workflow advisor by right-clicking on the *HDLTx*
% subsystem and selecting HDL Code -> HDL Workflow Advisor.
% # In step 1.1, at *Target workflow*, select *Software Defined Radio*. At
% *Target platform*, select *ZC706 and FMCOMMS1 RevB/C*. The
% *Synthesis tool* field should automatically set to *Xilinx Vivado* if you
% added Vivado to the path. Click apply.
% # In step 1.2, ensure that *User logic data path* is set to *Transmit
% path*. Click apply.
% # The default values for step 2 and step 3 can be used.
% # Run the Workflow Advisor to step 4. A quick way to run all the
% preceding steps is to right click on *4. Build SDR* in the workflow
% window and click *Run to selected task*.
%
% If no errors are found, the FPGA programming file generation process
% starts in an external command shell. You can monitor the implementation
% progress in the external shell. A message indicating successful
% completion of file generation is displayed in the shell upon completion.
%
% The HDL Workflow Advisor generates an SD card image in the folder
% _hdl_prj/sdr_prj/sdcard_image_, assuming the default folder was used in
% step 1.1. Use the following command to download the new image to the SD
% card:
%    
%    dev = sdrdev('ZC706 and FMCOMMS1 RevB/C');
%    downloadImage(dev, 'SDCardImage','C:\mywork\hdl_prj\sdr_prj\sdcard_image');
%
% The FPGA implementation of the QPSK transmitter can be verified using the
% companion retargeted model <matlab:sdrzQPSKTxFPGAretarget
% sdrzQPSKTxFPGAretarget>, whose transmit structure is shown in the figure
% below.
modelname = 'sdrzQPSKTxFPGAretarget';
load_system(modelname);
open_system(modelname);
%%
% This model can be built from the original <matlab:sdrzQPSKTxFPGA
% sdrzQPSKTxFPGA> model by removing the *Data Type Conversion* block and
% the *HDLTx* subsystem. *The functionality of these blocks is now
% implemented on the FPGA.*
%
% With the targeted SD card image loaded, simulating the
% sdrzQPSKTxFPGAretarget model will start transmitting the 'Hello world'
% messages.
close_system(modelname);


%% Transmitter Design: System Architecture
%
% The structure of the full QPSK transmitter in <matlab:sdrzQPSKTxFPGA
% sdrzQPSKTxFPGA> is shown in the figure below.
modelname = 'sdrzQPSKTxFPGA';
load_system(modelname);
open_system(modelname);
%%
% The HDL-optimized part of the QPSK transmitter is modeled in the *HDLTx*
% subsystem, whose contents is shown below. This is the DUT that will be
% implemented in the transmit path on the FPGA fabric.
open_system([modelname '/HDLTx']);
%%
% There are two implementation details that are specific to targeting the
% SDR hardware: the packing and unpacking of data at the input boundary of
% the *HDLTx* subsystem and the repeat/downsampling by 4.
% 
% * To conform to the SDR hardware interface, the input of the
% HDL-optimized QPSK transmitter must be a 16-bit signed complex value. To
% meet this requirement, the *Host to FPGA Packing* subsystem takes the
% incoming data bits and splits them into IQ channels. Each bit then
% becomes the LSB of a complex 16-bit integer, which is compatible with the
% SDR hardware interface. On the FPGA, the *Host to FPGA Unpacking*
% subsystem extracts the two LSBs and concatenates them into an unsigned
% 2-bit integer to be modulated by the *Symbol Mapping* subsystem.
% * The downsample by 4 at the input to the *HDLTx* subsystem is paired
% with the *Repeat* block at the top level of the model. When targeting a
% subsystem for inclusion in the SDR transmit chain, the subsystem must
% have the same sample rate at its input and output. Since the root raised
% cosine filter increases the sample rate by 4, it is necessary to
% artificially increase the sample rate at the input to the subsystem by 4.
% See <matlab:sdrzdoc('sdrz_targeting') Targeting Limitations>
% for more information.
close_system(modelname);


%% List of Example Helper Files
%
% This example uses the following helper files:
%
% * <matlab:edit('sdrzQPSKTxFPGA_init.m') sdrzQPSKTxFPGA_init.m>:
% initialize variables used for running the model.


displayEndOfDemoMessage(mfilename)