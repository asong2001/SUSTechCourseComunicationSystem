%% Targeting HDL Optimized QPSK Receiver Using Analog Devices AD9361/AD9364
%
% This example shows how to model an HDL-optimized QPSK receiver and
% prototype it on the FPGA of an SDR platform using the HDL Coder(TM) workflow
% advisor. The SDR device in this model is designed to receive indexed
% 'Hello world' messages and print them to the command window. For a full
% description of targeting, see <matlab:sdrzdoc('sdrz_targeting')
% Targeting Overview> in the documentation.
%
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting Started>
% documentation for details on configuring your host computer to work with
% the Support Package for Xilinx(R) Zynq-Based Radio. Additionally,
% targeting requires HDL Coder and Xilinx Vivado(R) 2015.2.1, which must
% be on the system path. You can add Vivado to the path when you call
% <matlab:sdrzdoc('sdrz_settoolpath') setupTools>. See the
% <matlab:sdrzdoc('sdrz_targeting') targeting requirements> for more
% details.

% Copyright 2014-2015 The MathWorks, Inc.

%% Introduction
%
% The <matlab:commqpskrxhdl commqpskrxhdl> model walks through building a
% practical HDL-optimized digital receiver, which includes coarse frequency
% compensation, PLL-based fine frequency compensation, timing recovery,
% frame synchronization, phase ambiguity resolution, and QPSK demodulation.
% This example takes the next steps to generate HDL code from this
% HDL-optimized receiver and prototype it on the FPGA of the SDR hardware.
%
% Targeting typically involves a pair of models:
%
% # The targeting model, used to generate the HDL code and integrate the
% DUT (Design Under Test) into the receive path implemented on the FPGA of
% the SDR hardware
% # The retargeted model, which uses the SDR hardware with the included DUT
% performing a portion of the baseband processing
%
% Note that there are two methods of running this example, based on the
% targeting/retargeted pair of models:
%
% # Run the <matlab:zynqRadioQPSKRxFPGAAD9361AD9364SL zynqRadioQPSKRxFPGAAD9361AD9364SL> targeting model
% directly. See the <#3 Running the Example: Simulated> section.
% # Generate a targeted bitstream and use the retargeted model
% <matlab:zynqRadioQPSKRxFPGAAD9361AD9364RetargetSL zynqRadioQPSKRxFPGAAD9361AD9364RetargetSL>. See the <#4
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
% 2. Ensure that you have a suitable transmitter. This example is designed
% to work in conjunction with any of the following possible transmitters:
%
% * The <matlab:showdemo('zynqRadioQPSKTransmitRepeatAD9361AD9364ML') QPSK Transmit Repeat Using Analog Devices AD9361/AD9364> MATLAB(R) example
% * The <matlab:showdemo('zynqRadioQPSKTxAD9361AD9364SL') QPSK Transmitter Using Analog Devices(TM) AD9361/AD9364>
% Simulink(R) example
% * The <matlab:showdemo('zynqRadioQPSKTransmitterAD9361AD9364ML') QPSK Transmitter Using Analog Devices AD9361/AD9364> MATLAB example
% * The <matlab:showdemo('zynqRadioQPSKTxFPGAAD9361AD9364SL') Targeting HDL-Optimized QPSK
% Transmitter Using Analog Devices AD9361/AD9364> Simulink example


%% Running the Example: Simulated
%
% To run the example directly, open <matlab:zynqRadioQPSKRxFPGAAD9361AD9364SL zynqRadioQPSKRxFPGAAD9361AD9364SL>
% and start the simulation. Note that this simulates the entire baseband
% design in a similar way to the <matlab:showdemo('zynqRadioQPSKRxAD9361AD9364SL') QPSK
% Receiver Using Analog Devices AD9361/AD9364> Simulink example i.e. no HDL code is
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
% The <matlab:zynqRadioQPSKRxFPGAAD9361AD9364SL zynqRadioQPSKRxFPGAAD9361AD9364SL> model can be used to generate
% a targeted SD card image for the SDR hardware that includes the *HDLRx*
% subsystem as a DUT. The DUT is automatically placed in the SDR receiver
% chain that is implemented on the FPGA.
%
% To generate the new SD card image, perform the following steps:
%
% # Make sure you have added Vivado to the system path. See
% <matlab:sdrzdoc('sdrz_settoolpath') setupTools>.
% # Open the HDL workflow advisor by right-clicking on the *HDLRx*
% subsystem and selecting HDL Code -> HDL Workflow Advisor.
% # In step 1.1, at *Target workflow*, select *Software Defined Radio*. At
% *Target platform*, select the appropriate option for your radio hardware (e.g. 'ZC706 and FMCOMMS2/3/4', 'Zedboard and FMCOMMS2/3/4' or 'PicoZed SDR'). The
% *Synthesis tool* field should automatically set to *Xilinx Vivado* if you
% added Vivado to the path. Click apply.
% # In step 1.2, ensure that *User logic data path* is set to *Receive
% path*. If you are using a Zedboard, you might need to set the synthesis frequency to a lower rate (1MHz would be appropriate). Click apply.
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
%    dev = sdrdev('ZC706 and FMCOMMS2/3/4');
%    downloadImage(dev, 'SDCardImage','C:\mywork\hdl_prj\sdr_prj\sdcard_image');
%
% The FPGA implementation of the QPSK receiver can be verified using the
% companion retargeted model <matlab:zynqRadioQPSKRxFPGAAD9361AD9364RetargetSL
% zynqRadioQPSKRxFPGAAD9361AD9364RetargetSL>, whose receiver structure is shown in the figure
% below.

modelname = 'zynqRadioQPSKRxFPGAAD9361AD9364RetargetSL';
load_system(modelname);
open_system([modelname '/QPSK Receiver']);
%%
% This model can be built from the original
% <matlab:zynqRadioQPSKRxFPGAAD9361AD9364SL
% zynqRadioQPSKRxFPGAAD9361AD9364SL> model by making the following
% modifications:
%
% * Removing the Data Type Conversion block, the Unbuffer block, the
% *HDLRx* subsystem, and the Buffer blocks. *The functionality of these
% blocks is now implemented on the FPGA.*
% * Deselecting the _Bypass user logic_ option in the _Advanced tab_ of the
% Zynq SDR Receiver block.
%
% With the targeted SD card image loaded, simulating the
% <matlab:zynqRadioQPSKRxFPGAAD9361AD9364RetargetSL
% zynqRadioQPSKRxFPGAAD9361AD9364RetargetSL> model should produce the
% expected 'Hello world' messages in the MATLAB command window.
close_system('zynqRadioQPSKRxFPGAAD9361AD9364RetargetSL/QPSK Receiver')
close_system('zynqRadioQPSKRxFPGAAD9361AD9364RetargetSL')


%% Receiver Design: System Architecture
% The structure of the QPSK receiver in <matlab:zynqRadioQPSKRxFPGAAD9361AD9364SL
% zynqRadioQPSKRxFPGAAD9361AD9364SL> is shown in the figure below.
modelname = 'zynqRadioQPSKRxFPGAAD9361AD9364SL';
load_system(modelname);
open_system([modelname '/QPSK Receiver']);


%%
% The HDL-optimized part of the QPSK receiver is modeled under the
% subsystem *HDLRx*, whose contents are shown below. This is the DUT that
% will be implemented in the receive path on the FPGA fabric.
close_system([modelname '/QPSK Receiver']);
open_system([modelname '/QPSK Receiver/HDLRx']);


%%
% Compared to the implementation of the *HDLRx* subsystem in the
% <matlab:commqpskrxhdl commqpskrxhdl> model, everything stays the same
% except for the presence of the *FPGA to Host* subsystem after *Data
% Decoding*. To conform with the Xilinx(R) FPGA interface in the SDR
% hardware, the outputs of the HDL-optimized QPSK receiver must be a 16-bit
% signed complex value. To meet this specific interface requirement, the
% subsystem *FPGA to Host* embeds the three Boolean signals (bit1, bit2,
% and dValid) into a 16-bit complex-valued integer.
%
% Corresponding to the presence of the *FPGA to Host* subsystem in *HDLRx*,
% the *Unpack FPGA Outputs* subsystem extracts the three Boolean signals
% (bit1, bit2, and dValid) from the 16-bit signed complex values. The
% extracted Boolean signals are then fed to the *dataframer* subsystem.
close_system([modelname '/QPSK Receiver/HDLRx']);
close_system(modelname);

%% List of Example Helper Files
%
% This example uses the following helper files:
%
% * <matlab:edit('zynqRadioQPSKRxFPGAAD9361AD9364SL_init.m') zynqRadioQPSKRxFPGAAD9361AD9364SL_init.m>:
% initialize variables used for running the model.


displayEndOfDemoMessage(mfilename)

