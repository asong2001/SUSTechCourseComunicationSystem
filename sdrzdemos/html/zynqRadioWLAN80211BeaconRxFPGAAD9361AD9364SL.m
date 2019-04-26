%% Targeting the HDL Optimized 802.11 WLAN Beacon Frame Receiver Using Analog Devices AD9361/AD9364
%
% This example shows how to run an HDL-optimized Beacon frame receiver
% model with  AD9361/AD9364 SDR hardware and how to prototype it on the
% FPGA of the SDR platform using the HDL Coder(TM) workflow advisor. For a
% full description of targeting, see <matlab:sdrzdoc('sdrz_targeting')
% Targeting Overview> in the documentation.
%
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting Started>
% documentation for details on configuring your host computer to work with
% the Support Package for Xilinx(R) Zynq-Based Radio. Additionally,
% targeting requires HDL Coder and Xilinx Vivado 2015.2.1, which must be on
% the system path. You can add Vivado to the path when you call
% <matlab:sdrzdoc('sdrz_settoolpath') setupTools>. See the
% <matlab:sdrzdoc('sdrz_targeting') targeting requirements> for more
% details.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
%
% The <matlab:showdemo('commwlan80211BeaconRx') IEEE 802.11 WLAN - Beacon
% Frame Receiver with Captured Data example> shows a model that receives
% beacon frames in an 802.11 wireless local area network (WLAN) as
% described in [ <#21 1>]. This example shows a hardware friendly version
% of the model and how to use it with SDR hardware, as well as the steps
% taken to prototype it on the FPGA.
%
% There are two methods of running this example, based on the
% targeting/retargeted pair of models:
%
% # Run the <matlab:zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL
% zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL> targeting model directly.
% See the <#14 Running the Example: Simulated> section.
% # Generate a targeted bitstream and use the retargeted model
% <matlab:zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364RetargetSL
% zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364RetargetSL>. See the <#16
% Running the Example: Targeted> section.
%
%% Zynq-Radio WLAN Beacon Receiver
% The hardware friendly version of the
% <matlab:showdemo('commwlan80211BeaconRx') IEEE 802.11 WLAN - Beacon Frame
% Receiver with Captured Data example> is placed within an enabled
% subsystem _802.11 WLAN Beacon Frame Receiver_, which ensures that only
% valid samples coming from the SDR hardware will be processed.

modelname = 'zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL';
open_system(modelname);
%% 
open_system([modelname '/802.11 WLAN Beacon Frame Receiver']);

%%
% The HDL-optimized subsystem is called *HDLRx* and this is the DUT that
% will be implemented in the FPGA fabric. The contents of the subsystem are
% displayed below.

open_system([modelname '/802.11 WLAN Beacon Frame Receiver/HDLRx']);
%%
%
% The *HDLRx* subsystem, in partnership with the *Unpack FPGA Outputs*
% subsystem of the model, ensures that the the start signal can be packed
% into the 16 bit complex value required by the hardware.
%
% The *FPGA to Host* subsystem packages the start flag and complex
% fixed-point data by replacing the least significant bit of the 'Symbols'
% output with the boolean value representing the start flag. The *Unpack
% FPGA Outputs* subsystem is then used to recover the start flag and
% complex data representing the recovered symbols for further processing in
% Simulink.
%% HDL Modifications to Beacon Frame Receiver Example
%
% The <matlab:showdemo('commwlan80211BeaconRx') IEEE 802.11 WLAN - Beacon
% Frame Receiver with Captured Data example> consists of three main
% components: front end, receiver controller, and detector. The front end
% and the receiver controller operate at a high rate in the receiver, so
% they have been optimized for HDL code generation in this example. The
% following sections describe the details of these modifications.
%
% To generate efficient HDL code the following modifications to the model
% are required:
% 
% * *Streaming Input and Output:* The HDL optimized beacon frame receiver
% processes data one sample at a time. The real-world signal is streamed
% into the receiver front-end. The streaming output of the receiver
% controller is buffered and passed to the detector, which operates on the
% data on a per-frame basis.
% 
% * *Fixed-point:* The receiver front-end and the controller logic operate
% in fixed-point mode.
% 
% * *HDL optimized architecture:* Several blocks have been redesigned to
% use hardware efficient algorithms and architectures.
%
modelname = 'zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL';

% Define Simulink(R) blocks as variables

agcScope = [modelname ...
    '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Front End/AGC Results/AGC Scope'];
freqOffsetScope = [modelname ...
    '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Front End/Frequency Offset'];
syncScope = [modelname ...
    '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Controller/Sync Results/Synchronization Scope'];
sfdSyncScope = [modelname ...
    '/802.11 WLAN Beacon Frame Receiver/Detector/SFD Synchronization Scope'];
symbolsSP = [modelname ...
    '/802.11 WLAN Beacon Frame Receiver/Frame Buffer/Received Symbols'];


%% Receiver Front-End
% The receiver front-end design is modeled in the subsystem *HDLRx Front
% End*. The receiver front-end is composed of a matched-filter, AGC, and
% coarse frequency compensation.

open_system([modelname '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Front End']);
%%
% Modifications were made to the coarse frequency estimation algorithm [
% <#21 2> ] implemented in the original beacon receiver model:
%
% * The auto correlation operation has been replaced by a simple smooth
% filter.
% * The _angle_ function has been implemented using the _Complex to
% Magnitude-Angle HDL Optimized_ block. This block computes the phase using
% the hardware friendly CORDIC algorithm. To learn more about the _Complex
% to Magnitude-Angle HDL Optimized_ block, refer to the DSP System Toolbox
% <matlab:helpview(fullfile(docroot,'toolbox','dsp','dsp.map'),'dsphdlmagangle')
% documentation>.
% * The detected phase offset is sent to _NCO HDL Optimized_ block to
% generate a complex exponential signal that is used to correct the phase
% offset in the original signal.
open_system([modelname '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Front End/Coarse Frequency Compensation/luise_fixpt'],'force');

%%
% The _NCO HDL Optimized_ block provides hardware friendly options, maps the
% lookup table into a ROM, and provides a lookup table compression option
% to significantly reduce the lookup table size. To learn more about HDL
% support for _HDL Optimized NCO_ block, refer to the
% <matlab:helpview(fullfile(docroot,'toolbox','dsp','dsp.map'),'nco_hdl')
% documentation>.
close_system([modelname '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Front End/Coarse Frequency Compensation/luise_fixpt']);
open_system([modelname '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Front End/Coarse Frequency Compensation']);

%% Receiver Controller 
% The receiver controller finds the correlation between the phase offset
% corrected signal and synchronization signal. When the peak of the
% correlation is detected, the signal is delayed based on the peak position
% before despreading.  Because of the large beacon frame size (2816
% samples), implementing the correlator using either an FFT or a matched
% filter is not efficient in hardware. In addition, finding the maximum of
% a 2816-element vector is not hardware friendly.  The correlation and
% despreading algorithm has been redesigned in this example.  Despreading
% is performed before the correlation, aiming to reduce the filter size
% from 2816 to 128.  Because the start of the beacon signal is unknown when
% performing the despreading early, 22 channels have been designed in the
% Despread_Matched Filter module, with each channel offsetting the adjacent
% channel by one sample. The maximum of the outputs of 22 filters is
% computed and the despread results from the channel that produces the
% maximum are selected to send to the Detector.
%

close_system([modelname '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Front End/Coarse Frequency Compensation/luise_fixpt']);
open_system([modelname '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Controller']);

%%
% The downsample operation before the 22-channel matched filter in the
% _Despread_Matched Filter_ subsystem shown below allows the possibility of
% sharing resources across the 22 channels of the FIR matched filters, as a
% way to balance the hardware speed and resource usage. To enable the
% sharing of resources, the _ChannelSharing_ HDL property for the Discrete
% FIR Filter has been turned on. this reuses the FIR Filter hardware across
% all 22 channels. To see this property, right-click on the Discrete FIR
% Filter block, choose _HDL Code_, then _HDL Block Properties_. From the
% command line you can use hdlget_param to get, hdlset_param to set and
% hdldispblkparams to display the HDL properties of a block.

close_system([modelname '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Controller']);
open_system([modelname '/802.11 WLAN Beacon Frame Receiver/HDLRx/HDLRx Controller/Despread_Matched Filter'],'force');

%%
% The despreaded signal is buffered in frames of 128 symbols before being
% processed in the detector.  In the original Beacon Frame Receiver
% example, the PLCP information was fed back to the receiver controller to
% decide the length of payload to be collected.  In this example, in order
% to simplify the interface between the hardware and software components
% and not require real-time communication of information between the two,
% the maximal length of the payload is collected in the hardware and sent
% to the software.


%% Setup
%
% Before running the example, ensure you have performed the following
% steps:
%
% 1. Configure your host computer to work with the Support Package for
% Xilinx Zynq-Based Radio. See <matlab:sdrzdoc('sdrzspsetup') Getting
% Started> for help.
%
% 2. Make sure you have a suitable source to receive the 802.11 Beacons
% from. This can be done in multiple ways:
%
% * Grab the signals off the air from a local access point. You will need
% to find out which channel it will be transmitting from. The most common
% channels in use are 1, 6 and 11.
% * Use a wireless router and optionally set the beacon interval to a
% convenient value (usually 50 or 100ms). Having a lower beacon interval
% time will allow more beacons to be captured in a burst and help configure
% the model.


%% Running the Example: Simulated
%
% To run the example directly, open
% <matlab:zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL
% zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL>, configure the model
% parameters and start the simulation.
%
% The *Model Parameters* block can be used to modify:
%
% * Radio parameters - the channel number, RF gain and center frequency
% offset.
% * Receiver parameters - AGC max gain and step size, correlation
% threshold. (These are used to tune the HDL-optimized receiver subsystem
% and will affect the hardware implementation).
%

%%
% Note that since the model is optimized for hardware, it will run
% significantly slower than its floating point, frame-based counterparts.
% Due to this, the example uses <matlab:sdrzdoc('sdrz_burstmode') burst
% mode> to receive chunks of contiguous data (otherwise the SDR would be
% passing samples through at a much faster rate than the Simulink model can
% process them, resulting in lost samples).
%
% The SDR hardware is set to capture 103ms worth of data per burst (2816
% samples per frame x 805 frames per burst / 22Msps). This should be enough
% to receive at least one beacon frame on each burst in a standard setting
% (most access points will be transmitting a beacon every 100ms with the
% beacon lasting 3ms, however beacon intervals can vary).
%
% To obtain the scope outputs seen below, a wireless router was configured
% to transmit with a beacon interval of 50ms.
%
% <<zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL_Scopes.png>>
%
% As seen from the synchronization scope, there were 2 distinct instances
% 50ms apart where the correlation threshold was reached. The SFD
% Synchronization scope outputs are only on for the duration of packet
% detection which is turned on by the 'start' flag from the *HDLRx
% Controller*. To read more about the scope outputs, visit the
% <matlab:showdemo('commwlan80211BeaconRx') IEEE 802.11 WLAN - Beacon Frame
% Receiver with Captured Data example>.
%
% The decoded packet service type can be seen in the *MAC Header*
% information tab of the MPDU GUI as shown below.
%
% <<zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL_Header.png>>
%
% The receiver also decodes the SSID, seen in the *Information Elements*
% tab.
%
% <<zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL_Info.png>>
%

%% Running the Example: Targeted
%
% The <matlab:zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL
% zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL> model can be used to
% generate a targeted SD card image for the SDR hardware that includes the
% *HDLRx* subsystem as a DUT. The DUT is automatically placed in the SDR
% receiver chain that is implemented on the FPGA.

%%
% To generate the new SD card image, perform the following steps:
%
% # Make sure you have added Vivado to the system path. See
% <matlab:sdrzdoc('sdrz_settoolpath') setupTools>.
% # Open the HDL workflow advisor by right-clicking on the *HDLRx*
% subsystem and selecting HDL Code -> HDL Workflow Advisor.
% # In step 1.1, at *Target workflow*, select *Software Defined Radio*. At
% *Target platform*, select the appropriate option for your radio hardware
% (e.g. 'ZC706 and FMCOMMS2/3/4', 'Zedboard and FMCOMMS2/3/4' or 'PicoZed
% SDR'). The *Synthesis tool* field should automatically set to *Xilinx
% Vivado* if you added Vivado to the path. Click apply.
% # In step 1.2, ensure that *User logic data path* is set to *Receive
% path*. If you are using a Zedboard, you might need to set the synthesis
% frequency to a lower rate (30MHz would be appropriate). Click apply.
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
% hdl_prj/sdr_prj/sdcard_image/, assuming the default folder was used in
% step 1.1. Use the following command to download the new image to the SD
% card:
%    
%    dev = sdrdev('ZC706 and FMCOMMS2/3/4'); downloadImage(dev,
%    'SDCardImage','C:\mywork\hdl_prj\sdr_prj\sdcard_image');
%
%%
%
% The FPGA implementation of the beacon receiver can be verified using the
% companion retargeted model
% <matlab:zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364RetargetSL
% zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364RetargetSL>.
%
% This model can be built from the original
% <matlab:zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL
% zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL> model by removing the *Data
% Type Conversion* and *HDLRx* blocks from the subsystem. *The
% functionality of these blocks will be implemented on the FPGA.*

modelname = 'zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364RetargetSL';
load_system(modelname);
open_system([modelname '/802.11 WLAN Beacon Frame Receiver']);
%%
%
% Running this model, with the SD Card image downloaded, you should notice
% the model performing much more quickly, as now most of the baseband
% processing is being performed on the FPGA. Also, since the HDL-optimized
% subsystem is sending back downsampled, more tightly packed information in
% the same data burst, there will be many more beacons demodulated and
% decoded per burst as shown by the SFD Synchronization scope below.
%
% <<zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364RetargetSL_SFD.png>>
%
% *Note that the scope no longer represents 103ms of recorded data from the
% RF card. The new time interval of captured and processed data will be 22
% times longer at 2266.9ms (2816 samples per frame x 805 frames per burst /
% 1Msps).*
%
% Another important fact to consider is that the *HDLRx* subsystem
% parameters for AGC and Synchronization will now be hard coded on the FPGA
% fabric, thus it is advisable to first run and tune the simulation model
% before attempting to target the HDL-optimized subsystem.
%


%% List of Example Helper Files
%
% This example uses the following helper files:
%
% * <matlab:edit('zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL_params.m')
% zynqRadioWLAN80211BeaconRxFPGAAD9361AD9364SL_params.m>: initialize
% variables used for running the model.


%% References
% # IEEE Std 802.11-2007: _IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements, Part 11: Wireless LAN
% Medium Access Control (MAC) and Physical Layer (PHY) Specifications_, 
% IEEE, New York, NY, USA, 1999-2007.
% # M. Luise and R. Reggiannini, "Carrier frequency recovery in all-digital
% modems for burst-mode transmissions," _IEEE Trans. Communications_, pp.
% 1169-1178, Feb.-March-Apr. 1995.

displayEndOfDemoMessage(mfilename)

