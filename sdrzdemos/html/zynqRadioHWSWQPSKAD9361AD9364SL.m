%%  HW/SW Co-design QPSK Transmit and Receive Using Analog Devices AD9361/AD9364 
%
% This example shows how to implement algorithms on the Zynq radio platform
% that are partitioned across the ARM and the FPGA fabric. QPSK transmit
% and receive functions are implemented in hardware and software and mapped
% to the radio platform as shown in the diagram below.
%
% <<zynqRadioHWSWQPSKAD9361AD9364SL_QPSKBlockDiagram.png>>
%
% Required products:
%
% * Simulink
% * Communications System Toolbox
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
% This example is based on existing
% <matlab:showdemo('zynqRadioQPSKTxFPGAAD9361AD9364SL') transmit> and
% <matlab:showdemo('zynqRadioQPSKRxFPGAAD9361AD9364SL') receive> QPSK
% targeting examples present in this support package. The transmit and
% receive FPGA implementations are combined into one HDL IP core and
% implemented on the Zynq programmable logic (PL). The data encoding and
% decoding, which in the other examples is undertaken on the host machine,
% is here run on the Zynq ARM processor through code generation. Some
% control parameters are added to the FPGA IP core to show how the design
% can be adjusted in real time using AXI4-Lite registers accessed from
% Simulink.

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
%
% The hardware generation model is used to develop the functionality that
% you wish to implement on the FPGA fabric. This is similar to the path
% taken in the <matlab:showdemo('zynqRadioQPSKRxFPGAAD9361AD9364SL')
% Targeting HDL Optimized QPSK Receiver Using Analog Devices AD9361/AD9364>
% and <matlab:showdemo('zynqRadioQPSKTxFPGAAD9361AD9364SL') Targeting HDL
% Optimized QPSK Transmitter Using Analog Devices AD9361/AD9364> examples,
% where an HDL-optimized QPSK transmitter and receiver are modelled and
% then implemented on the FPGA fabric of the Zynq radio platform. In this
% example, we transmit and receive on a single board, however, the example
% could be modified to work in frequency division duplex (FDD) by moving
% the transmit and receive center frequencies apart and using the same
% model on two separate boards.
%
% <matlab:open_system('zynqRadioHWSWQPSKAD9361AD9364SL') Open the model.>
% 
open_system('zynqRadioHWSWQPSKAD9361AD9364SL')
%%
% *Hardware/software partitioning:* In general, the programmable logic of
% the FPGA is used for high rate signal processing while the ARM is used
% for slower rate, control functionality. In this example, the QPSK
% transmit and receive physical layer front ends are implemented on the
% programmable logic, as these include high rate operations such as gain
% control, filtering and frequency compensation. The data encoding and
% decoding are much slower rate and are implemented on the ARM, which
% outputs the decoded message to a console window.
% 
% The functionality used in this example is taken from the existing
% transmitter and receiver targeting examples, with some modifications. The
% IP Core Generation Workflow used to implement the FPGA IP core and
% generate the software interface model has some specific requirements:
% 
% * Complex inputs and outputs are not supported at the ports of the HDL
% subsystem. Therefore, real and imaginary signals must be modelled at the
% subsystem boundaries.
% * The data inputs and outputs to the subsystem are modelled using
% separate data and valid signals. The input and output clock rates of the
% subsystem must be equal, in Simulink the data and valid lines must be
% driven at the same sample rate. The valid signal must also be modelled at
% the input and output of the user logic.
%
% *AXI4-lite Control Port txSrcSelect:* A control input |txSrcSelect| is
% added to allow some control over the transmitted data. The |txSrcSelect|
% port on the HDL subsystem is used to select between two different data
% sources for the transmitter. If the |txSrcSelect| port is true, the data
% source for the transmitter will be a look-up table stored on the FPGA
% fabric and the received data should resemble the "Hello World 0XX"
% strings seen in the other QPSK examples. If the |txSrcSelect| port is
% false, the data source for the received data will be the ARM processor,
% which will generate samples in real-time and send them to the transmitter
% on the FPGA fabric. The message in this case will be "*Zynq HW/SW
% Co-design*". The message text is taken from a workspace variable |txtStr|
% which can be modified at compile time to change the message. Note that
% the length of this string must be less than 24 characters.
%
% You can run this model and confirm its operation. By double clicking the
% switch *Internal Tx Switch* you can select which source to transmit from.
% Note that *Goto* and *From* blocks are used to model the antenna
% connection and pass the transmitted data at the output of the transmit
% user logic to the input of the receive user logic.

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
% Start the targeting workflow by right clicking the 
% |HDL_QPSK| subsystem and selecting |HDL Code / HDL Workflow Advisor|.
%
%
% * In Step 1.1, select |IP Core Generation| workflow and the appropriate
% Zynq radio platform from the choices: |PicoZed SDR|, |ZC706 and
% FMCOMMS2/3/4|, |ZedBoard and FMCOMMS2/3/4|.
%
% <<zynqRadioHWSWQPSKAD9361AD9364SL_step1dot1.png>>
%
% * In Step 1.2, select |Receive and transmit path| reference design. The
% interface table can then be used to map the user logic signals to the
% interface signals available in the reference design. In this example, we
% are only using a single channel, so the channel 1 connections should be
% connected to the relevent ports as shown below. Ensure that the Channel
% Mapping parameter is set to 1, and that the DUT Synthesis Frequency is
% set to a reasonable number given the baseband sampling rate of the
% system. In the shipping example, the sample rate is just above 520ksps,
% so a synthesis frequency of 1MHz is sufficient.
%
% <<zynqRadioHWSWQPSKAD9361AD9364SL_step1dot2.png>>
%
%
close_system('zynqRadioHWSWQPSKAD9361AD9364SL',0)
%%
% * Step 2 prepares the model for HDL code generation by doing some design
% checks.
% * Step 3 performs the actual HDL code generation for the IP core.
% * Step 4 integrates the newly generated IP core into the larger Zynq SDR
% reference design, generates the bitstream and helps you load it onto the
% board.
%
% Execute each step in sequence to experience the full workflow, or, if you
% are already familiar with preparation and HDL code generation phases,
% right click Step 4.1 in the table of contents on the left hand side and
% select |Run to selected task|. You should not have to modify any of the
% default settings in Steps 2 or 3.
%
%% Software generation model and block library
% * In Step 4.2, the workflow generates a Zynq software generation
% interface model and a block library. Click the |Run this task| button
% with the default settings.
%
% <<zynqRadioHWSWQPSKAD9361AD9364SL_step4dot2.png>> 
% 
% *Software Interface Library*
% 
% <<zynqRadioHWSWQPSKAD9361AD9364SL_GenLibrary.png>>
% 
% The library contains the AXI Interface block which has been generated
% from the HDL_QPSK subsystem. Note that this exposes only the AXI4-lite
% control ports, and not the data ports. The data ports are present on the
% Receiver/Transmitter blocks which represent the data interface between
% the FPGA user logic and the ARM. If you use the library blocks in your
% downstream models, any updates you make to your HDL subsystem will
% automatically be propagated to this library and then to your software
% generation models when you run through the workflow. In this example, the
% hardware generation model did not contain any SDR transmit or receive
% blocks so the parameters on these blocks could not be populated. When
% using the library blocks you must ensure to configure the parameters
% correctly for your application.
% 
% *Software Interface Model*
% 
% <<zynqRadioHWSWQPSKAD9361AD9364SL_GenTemplate.png>>
% 
% The software interface model can be used as a starting point for full
% SW targeting to the Zynq: External mode simulation, Processor-in-the-loop
% and full deployment. Note that this generated model will be overwritten
% each time Step 4.2 is run, so it is advisable to save this model under
% a unique name and develop your software algorithm in there. A software
% interface model has been provided which shows how you may decide to
% structure this model, see section *Running the Software and Hardware on
% the Zynq board*.
%
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

%% Constructing the QPSK software interface model
% A software interface model has been provided which shows how you could
% modify the generated model to set it up for the QPSK example. This
% interface model will allow you to run the model in |External mode| or
% fully deployed.
%
% <matlab:open_system('zynqRadioHWSWQPSKAD9361AD9364SL_interface') Open the model.>

open_system('zynqRadioHWSWQPSKAD9361AD9364SL_interface')

%% Setting up the software model to run on the ARM processor
% The application model has been set up following the guidelines in the
% <matlab:sdrzdoc('hwswcodesign_workflow') HW/SW co-design workflow
% documentation>, section *Configure Software Interface Model*.
% 
% * The model is continuously transmitting and receiving data, so it has
% been configured to run from the Transmit interrupt. This ensures that the
% ARM and the FPGA are running in synchronisation and means that the
% software will be driven by a schedule tick at the frame rate.
%
% <<zynqRadioHWSWQPSKAD9361AD9364SL_ConfigParams.png>>
%
% * Buffer data for continuous transmission is selected on the transmit
% block, this ensures two frames of data are buffered before the
% transmitter is started. To work in this mode, the transmit frame size has
% been increased to 10000 by concatenating 50 frames at a time using a For
% Iterator and MATLAB Function block.
% * The transmitter underflow has not been connected to a Stop block as
% underflows may happen when the source switch is toggled when running in
% |External Mode|. This is as a result of the increased processing required
% by |External mode| register writes. See the section below on running the
% model using UDP blocks to control the hardware.
% * The receiver is placed within an enabled subsystem to delay it starting
% by 2 frame periods. This ensures the transmitter is running before the
% receiver for most robust performance.
% * As the QPSK receiver contains a downsample by two operation, the valid
% signal at the output of the user logic is used to reduce the effective
% sample rate to half that of the clock rate. In the model, the frame rate
% of the ARM processor is therefore set to half the sampling rate
% multiplied by the frame size. The frame size for the receiver is half
% that of the transmitter for this same reason.
% * The Receive Timeout has been set to 200 $\mu$ s, which is small in
% relation to the frame period.
%% Running the Software and Hardware on the Zynq board
%
% |External mode| allows you to control the configuration from the Simulink
% model. Once the design is running, switch between sourcing data from the
% ARM or the FPGA fabric by toggling the Tx Switch. 
% 
% You can also fully deploy the design to run on the board, disconnected
% from Simulink.  In the Simulink toolbar, click *Deploy to Hardware*. In
% this mode you will not be able to tune parameters.
%
close_system('zynqRadioHWSWQPSKAD9361AD9364SL_interface',0)

%% Using UDP blocks to control the user logic
% |External mode| requires some overhead to be included in the software
% running on the hardware to deal with communication between the host and
% the board. As seen in the software interface model, switching between the
% ARM and FPGA transmitter sources resulted in transmit underflows. An
% alternative interface model has been supplied which shows how UDP blocks
% can be used as an alternative switching mechanism which requires less
% overhead.
%
% <matlab:open_system('zynqRadioHWSWQPSKAD9361AD9364SL_UDP_interface') Open the model.>

open_system('zynqRadioHWSWQPSKAD9361AD9364SL_UDP_interface')

%%
% In this model, the switch has been replaced by a UDP receive block which
% will be able to receive UDP packets and output the source choice value.
% Some further modifications have been made to the model.
%
% * The transmitter underflow has now been connected to a Stop block which
% will cause the model to exit whenever an underflow is detected. A Step
% source has been used to gate the overflow signal for the first ten frame
% periods as the hardware starts up. Some underflows may be experienced as
% the transmitter performs buffering, and then again when the receiver
% starts two frames later. The receiver initialisation will result in a
% one-time load on the processor that may cause underflows.
% * The UDP Receive block has a default output value of 0, so an inverter
% has been placed between the UDP Receive block and the multiplexors to
% ensure that the transmitter starts up with the FPGA transmit source. In
% the Simulink toolbar, click *Deploy to Hardware*. Once the model has
% deployed to the hardware, a UDP transmitter source can be used to drive
% the transmit source selection.
close_system('zynqRadioHWSWQPSKAD9361AD9364SL_UDP_interface',0)
%%
% A simple UDP transmitter model has been supplied which can be used to
% drive the transmit source select.
%
% <matlab:open_system('zynqRadioHWSWQPSKAD9361AD9364SL_UDP_SourceSelect') Open the model.>

open_system('zynqRadioHWSWQPSKAD9361AD9364SL_UDP_SourceSelect')

%%
% This model has been configured to run for a single step and send a single
% UDP packet containing the source select value. Set the source select to
% your desired source and click play to send a UDP control packet.
close_system('zynqRadioHWSWQPSKAD9361AD9364SL_UDP_SourceSelect',0)

