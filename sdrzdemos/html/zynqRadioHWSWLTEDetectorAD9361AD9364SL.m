%% HW/SW Co-design Implementation of LTE OFDM Detector Using Analog Devices AD9361/AD9364
%
% This example shows how to implement an LTE OFDM detector on the Zynq radio
% platform that is partitioned across the ARM and the programmable logic
% (PL) fabric.
% 
% Required products:
%
% * Simulink
% * Communications System Toolbox
% * LTE System Toolbox
% * HDL Coder
% * HDL Coder Support Package for Xilinx Zynq-7000 Platform
% * Embedded Coder
% * Simulink Coder
% * Embedded Coder Support Package for Xilinx Zynq-7000 Platform
% * Communications System Toolbox Support Package for Xilinx Zynq-Based
% Radio (this package)
% 
%%
%
% Copyright 2016 The MathWorks, Inc.

%% Introduction
%
% The <matlab:hdlcoder_lteofdm_modDetect hdlcoder_lteofdm_modDetect> model
% is a hardware friendly model that is capable of detecting the cell
% identity of an LTE OFDM signal. This example shows how to modify the
% model to use it with SDR hardware, and the steps taken to prototype it on
% the Zynq platform.

%% Setup
%
% If you have not yet done so, run through the Target Updater portion of
% the Zynq Embedded Coder support package and the Zynq Radio support
% package installations by following the
% <matlab:sdrzdoc('hwswcodesign_install') instructions in the
% documentation>. You can launch the setup wizards again without running
% the full installation by calling the targetupdater command.
% 
%   >> targetupdater
%
% Setup the Zynq SDR reference designs by calling the setup function. This
% must be called once per MATLAB session before you first open Simulink.
%
%   >> setupzynqradioipcoregen
%
%% Hardware Generation Model
% The hardware generation model is used to develop the functionality that
% you wish to implement on the PL fabric. In this example we implement an
% LTE OFDM Detector based on the
% <matlab:showdemo('hdlcoder_lteofdm_modDetect') HDL Implementation of LTE
% OFDM Modulator and Detector> example.
%
% <matlab:open_system('zynqRadioHWSWLTEDetectorAD9361AD9364SL') Open the model>.

modelname = 'zynqRadioHWSWLTEDetectorAD9361AD9364SL';
load_system(modelname);
open_system(modelname);

%%
% *Hardware/software partitioning:* In this example, the detector algorithm
% is implemented entirely in the PL, due to the high rate signal processing
% requirements of the design. These include filtering, and FFT operations.
% The ARM is used for post detection calculation and output of the detected
% LTE cell identity.
%
% The modified LTE OFDM Detector is shown below.

currentSubSys = 'HDL_LTE';
open_system([modelname '/' currentSubSys]);

%%
% This is the same HDLRx subsystem as seen in the
% <matlab:hdlcoder_lteofdm_modDetect hdlcoder_lteofdm_modDetect> model with
% the addition of the *FPGA_to_Host* subsystem, the *FIR Decimation* block
% and some further complex signal routing and rate change operations.
%
% The complex-valued input and output signals of the detector subsystem have
% been split into separate real and imaginary signals. This is a
% requirement of the IP Core Generation Workflow that is used to generate
% the FPGA IP core. 
%
% In order to take advantage of hardware resource sharing to implement some
% of the filter architectures found in the detector, the data rate into the
% system is higher than required. The input data rate is 61.44 MHz whereas
% the required data rate is 1.92 MHz, which translates to an overclocking
% factor of $N = 32$. Therefore, upon entering the detector, the data is
% immediately downsampled by a factor of 32. An *FIR Decimator* block is
% used to implement a lowpass decimation filter which captures the central
% part of the received LTE waveform and downsamples the data to 1.92 MHz.
%
% As the FIR filters in the *PSS_Detection* subsystem are implemented using
% a partly-serial hardware architecture, hardware resource sharing is used
% to clock them at a higher clock rate than the required data rate. This
% allows the filters to be implemented using fewer resources than if they
% were implemented using a fully-parallel hardware architecture.
%
% The PSS matched filters are implemented as discrete FIR filters with a
% partly-serial architecture. The serial partition is specified as [32 32 32 32],
% i.e. 4 partions of length 32. In the FPGA implementation, this requires a
% clock rate of 61.44 MHz, i.e. 32 times faster than the filter sample
% rate of 1.92 Msps. For more information on configuring partly serial
% filters, see the documentation for
% <matlab:web(fullfile(docroot,'hdlcoder/ug/configuring-hdl-filter-architectures.html'))
% Configuring HDL Filter Architectures>.
%
% At the output, the *FPGA_to_Host* subsystem packages the LTE cell ID data
% into the first 9 bits of the real-valued signal component. The
% imaginary-valued component comprises of the primary cell ID in the lowest
% 2 bits, and the secondary cell ID in the next 8 bits. A multiplexed
% signal is then created by alternatively outputing samples from the
% packaged Cell ID signal and the detected frequency offset signal. The
% data is later unpacked, decoded and displayed by software running on the
% ARM.
%
% Another requirement of the IP Core Generation Workflow is that the input
% and output rates of the HDL subsystem must be the same. Therefore the
% data rate at the output is increased to 61.44 MHz using *Repeat* blocks.
% In order to prevent overloading the ARM with a continuous stream of data
% at 61.44 MHz, the _validOut_ signal is only driven high for 5000 samples
% every second. This allows the ARM to receive a frame of 5000 samples from
% the HDL subsystem once per second. The control of the _validOut_ signal
% is contained in the *Divide_Valid* subsystem.


%%
% You can run this model to confirm its operation using the generated LTE
% waveforms in *zynqRadioLTETransmitData.mat*. The MAT-file contains four
% LTE waveforms, each generated with a different Cell ID. The Cell ID
% specific to each waveform is provided in the variable name, i.e.
%  |zynqRadioLTETransmitData_CellID_16| was generated using Cell ID *16*,
% and so on. 
%
% As the model contains a large number of HDL-optimized blocks, requiring
% simulation using sample-based signals, the simulation runs very slowly.
% In order to verify that the correct LTE cell ID is being detected, the
% _validOut_ signal is not being used to enable the *Decode and Print*
% software algorithm and, instead, a the enable signal is being held
% constantly high. *NOTE:* this is for simulation in Simulink only. When
% running with the detector targeted to the Zynq platform, the _validOut_
% signal should be used to enable the software algorithm.
%

%% IP Core Generation Workflow
% Once you are satisfied with the simulation behaviour of the hardware
% subsystem, you can start the process of generating the HDL IP Core,
% integrating it with the SDR reference design and generating software to
% run on the ARM.
%
% In preparation for targeting, you must set up the Xilinx tool chain by
% invoking |hdlsetuptoolpath|.  For example:
%
%   >> hdlsetuptoolpath('ToolName','Xilinx Vivado','ToolPath','C:\Vivado\2015.2\bin\vivado.bat');
%
% Start the targeting workflow by right-clicking the *HDL_LTE* subsystem
% and selecting |HDL Code / HDL Workflow Advisor|.
%
% * In Step 1.1, select |IP Core Generation| workflow and the appropriate
% Zynq radio platform from the choices: |PicoZed SDR|, |ZC706 and
% FMCOMMS2/3/4|, |ZedBoard and FMCOMMS2/3/4|.
% 
% * In Step 1.2, select |Receive path| reference design. The interface
% table can then be used to map the DUT signals to the interface signals
% available in the reference design. In this example, we are only using a
% single channel, so the channel 1 connections should be connected to the
% relevent ports as shown below. Ensure that the Channel Mapping parameter
% is set to 1, and that the DUT Synthesis Frequency is set to a reasonable
% number given the baseband sampling rate of the system. In this example
% the sample rate is 61.44 MHz, so a synthesis frequency of 61.44 MHz is
% required.
%
% <<zynqRadioHWSWLTEDetectorAD9361AD9364SL_HDLWAStep1_2.png>>
%


%%
% * Step 2 prepares the design for generation by doing some design checks.
% * Step 3 performs the actual HDL code generation for the IP core.
% * Step 4 integrates the newly generated IP core into the larger Zynq SDR
% reference design, generates the bitstream and helps you load it onto the
% board.
%
% In Step 4.1, select |Compile Optimized| as the _Synthesis objective_.
% This tells Vivado which synthesis strategy to use. By relaxing
% the synthesis effort, in this case by selecting |Compile Optimized|, the
% design will synthesize succesfully.
%
% <<zynqRadioHWSWLTEDetectorAD9361AD9364SL_HDLWAStep4_1.png>>
%
% Execute each step in sequence to experience the full workflow, or, if you
% are already familiar with preparation and HDL code generation phases,
% right click Step 4.1 in the table of contents on the left hand side and
% select |Run to selected task|. You should not have to modify any of the
% default settings in Steps 2 or 3.

%% Bitstream generation and loading
% The rest of the workflow is used to generate a bitstream for the FPGA
% fabric and download it to the board.
%
% * In Step 4.3, the workflow advisor generates a bitstream for the FPGA
% fabric. You can choose to execute this step in an external shell by
% ticking the selection |Run build process externally|.  This selection
% allows you to continue using MATLAB while the FPGA image is being built.
% The step will complete in a couple of minutes after some basic project
% checks have been completed, and the step will be marked with a green
% checkmark. However, you must wait until the external shell shows a
% successful bitstream build before moving on to the next step.
%
% * In Step 4.4, you can download the completed bitstream to the target and
% reboot the board. If you are downloading over Ethernet, make sure you have
% configured the board IP address first.
% 
%   >> devzynq = zynq(); 
%   >> devzynq.IPAddress = '192.168.3.2'; % This address should match that of the board 
%
% You can also do this at the command line as follows.
% 
%   >> devzynq = zynq(); 
%   >> devzynq.IPAddress = '192.168.3.2'; % This address should match that of the board 
%   >> devzynq.programFPGA('hdl_prj\vivado_ip_prj\vivado_prj.runs\impl_1\system_wrapper.bit') % Path to the generated bitstream
%
% Note that this FPGA image will not be persistent across power cycles. If
% you want the generated FPGA image to be loaded each time you turn on the
% board, call the device method |downloadImage|:
% 
%   >> dev = sdrdev('PicoZed SDR');
%   >> downloadImage(dev,'FPGAImage','hdl_prj\vivado_ip_prj\vivado_prj.runs\impl_1\system_wrapper.bit');
% 
% This call renames the generated system_wrapper.bit to system.bin and
% downloads the file over an Ethernet connection to the radio hardware.
%% Running the Software and Hardware on the Zynq board
%
% A software interface model has been provided which will allow you to run
% the example system in |External mode| or fully deployed. You must select
% the correct SDR receiver block for the hardware you are using.
%
% <matlab:open_system('zynqRadioHWSWLTEDetectorAD9361AD9364SL_interface') Open the model>.

close_system(modelname);
modelname = 'zynqRadioHWSWLTEDetectorAD9361AD9364SL_interface';
load_system(modelname);
open_system(modelname);

%%
% 
% The software interface model has been configured to run from the receive
% interrupt as the data coming from the FPGA is all the software needs to
% deal with; there are no AXI reads or writes.
%
% <<zynqRadioHWSWLTEDetectorAD9361AD9364SL_ConfigParams.png>>
%
% In order to successfully detect an LTE cell ID with this example, you can
% run in one of two ways:
% 
% # Detect a live signal off-the-air. For this option you will need to know
% the transmission center frequency of and LTE cell tower in your area. The
% default center frequency for this example is 816 MHz.
% # Use another radio source for the transmission of a generated LTE
% waveform. You could use the waveforms in *zynqRadioLTETransmitData.mat*
% as the transmission source. 
% 
% To run from saved data on another SDR board, use the following commands
% at the command line to send the test data.
%
%   >> load zynqRadioLTETransmitData.mat
%   >> tx = sdrtx('PicoZed SDR','BasebandSampleRate',61.44e6);
%   >> tx.CenterFrequency = 816e6; % Set center frequency to same as receiver
%   >> transmitRepeat(tx,zynqRadioLTETransmitData_CellID_123); 
%
% |External mode| allows you to control the configuration from the Simulink
% model. You can also fully deploy the design to run on the board, disconnected
% from Simulink.  In the Simulink toolbar, click *Deploy to Hardware*.
%
% The HDL user logic will only return regular valid signals to the ARM once
% a valid signal has been received. As the receiver is running from the
% receive interrupt, it will only print information if this has happened.
% In the lack of a valid LTE signal, this can also result in the ARM
% software locking up as it is not receiving any schedule ticks. You can
% release the software on the ARM using the following commands.
%
%   >> devzynq = zynq(); 
%   >> devzynq.IPAddress = '192.168.3.2'; % This address should match that of the board 
%   >> devzynq.stop('zynqRadioHWSWLTEDetectorAD9361AD9364SL');
%
% You could also consider configuring the model to run on a timer schedule
% which would allow it to keep running in the absence of data.

%%
% When running successfully, you should notice that the information printed
% to the console window updates once per second as expected. An example of
% the console output is shown below, where a Cell ID of *401* has been
% detected. The primary cell identity (_PCell ID_) and secondary cell
% identity (_SCellID_), along with the real-time  carrier frequency offset
% (_Freq Est_) are also shown. The hardware user logic works on a one shot
% basis, so tuning the centre frequency in real time will not result in
% detection of new LTE Cell IDs.
%
% <<zynqRadioHWSWLTEDetectorAD9361AD9364SL_LiveCapture.png>>

close_system(modelname);
displayEndOfDemoMessage(mfilename)
