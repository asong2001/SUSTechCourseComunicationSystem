%% IEEE 802.11b WLAN - Beacon Frame Receiver Using Analog Devices AD9361/AD9364
%
% This example shows reception of beacon frames in an 802.11b wireless
% local area network (WLAN) as described in [ <#11 1> ]. The example
% utilizes SDR hardware to receive radio signals and transfer them to
% Simulink(R) for processing. For more information refer to
% <matlab:showdemo('commwlan80211Beacon') IEEE(R) 802.11 WLAN - Beacon
% Frame> and <matlab:showdemo('commwlan80211BeaconRx') IEEE 802.11
% WLAN - Beacon Frame with Captured Data> examples.
%
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting
% Started> documentation for details on configuring your host computer to
% work with the Support Package for Xilinx(R) Zynq-Based Radio.
%
% Copyright 2014 The MathWorks, Inc.

%% Structure of the Example
% The model has three main parts: 
% 
% * |Model Parameters| block, where you can adjust several receiver
% parameters
% * 802.11b receiver, which comprises a receiver front end, receiver
% controller, and detector
% * Results, where you view several signals and the received information
%
% The following sections describe modifications made to the model presented
% in <matlab:showdemo('commwlan80211BeaconRx') IEEE 802.11 WLAN -
% Beacon Frame Receiver with Captured Data> example to make it work with
% the SDR hardware.
%%
open_system('zynqRadioWLAN80211BeaconRxAD9361AD9364SL');

%%
% This 802.11b WLAN example includes all the receiver signal processing in
% an enabled subsystem. Connecting the |DataLength| output of the SDR
% receiver block to the enable input of the subsystem ensures that the
% receiver only processes valid data.


%% Running the Example
% You can observe several signals in the scopes. The MPDU (MAC Protocol
% Data Unit) GUI figure shows the PLCP (Physical Layer Convergence
% Procedure) and MPDU CRC status and also the content of correctly decoded
% MPDU packets.
%
% Comment out the scatter plots in |Receiver/Receiver Controller| if not
% needed (right mouse click, Comment out) to increase the simulation speed.
%
% Several parameters that influence the ability to receive the beacon are
% in the |Model Parameters| block:
%
% * If the received signal is too weak or too strong, you might notice some
% garbled message output. In that case, you can change the gain of the SDR
% receiver block for better reception via |RF board gain|.
% * You can determine the responsiveness of the automatic gain control
% using |AGC loop gain| and |Maximum AGC gain|.
% * You can adjust for slight center frequency mismatches between the
% transmitter and receiver using |RF center frequency offset|.
% * You can increase the signal magnitude going into the downstream
% processing using |Receiver frontend gain|.
%
% It is also worth noting that this example pushes the gigabit Ethernet
% link to near its limit. If you run the example and the scopes don't show
% any data being received at all, try running a less demanding example
% first. This will help identify if the problem is performance related or a
% hardware setup issue. The <matlab:showdemo('zynqRadioQPSKTxAD9361AD9364SL') QPSK
% Transmitter> and <matlab:showdemo('zynqRadioQPSKRxAD9361AD9364SL') QPSK Receiver> examples
% are a good starting point.


%% SDR Receiver
% 802.11b uses 1e6 symbols per second for beacon signaling. Since the
% standard [ <#12 1> ] calls for a spreading factor of 11, the chip rate
% is 11e6 chips per second. The receiver needs at least two samples per
% chip, which is 22e6 samples per second. The |Baseband sample rate| in the
% SDR Receiver block can be set as 22e6 samples per second,
% which is the required rate.
%
%%
open_system('zynqRadioWLAN80211BeaconRxAD9361AD9364SL/Receiver');

%%
% Running this receiver simulation requires more time than processing the
% same data in real-time, especially when using the visualization scopes.
% To help alleviate this time requirement, the SDR receiver block uses
% <matlab:sdrzdoc('sdrz_burstmode') burst mode processing>.  Burst mode
% processing enables you to utilize the visualization capabilities of
% Simulink, while processing real data without the need of capturing and
% saving it.
%
% In burst mode, the block stores a contiguous burst of samples. The number
% of samples is determined by the values specified in the |Frame length| 
% and |Number of frames in burst| parameters. Each Simulink
% time step, the SDR receiver block sends a frame of samples to the
% Receiver subsystem. Most Wi-Fi(R) routers use a beacon interval of 100 Time
% Units (TU), which is 102.4 msec and the beacon packet lasts
% approximately 3 msec. Therefore, the receiver requires at least 106 msec
% of data to receive one beacon packet.


%% Exploring the Example
% You can try different channel numbers from the |Model Parameters| block
% mask. The most widely used channels are 6 and 11.
%
% This example allows you to modify several receiver parameters through the
% |Model Parameters| block mask dialog to optimize the receiver
% performance. If you notice that the AGC (Automatic Gain Control) gain
% reaches its maximum gain even when your signal is present at the receiver
% input, increase the maximum gain of the AGC. If the AGC is slow to
% respond to changes in the input signal amplitude, increase the AGC loop
% gain. Observe the AGC behavior in the AGC scope.
%
% If your signal results in smaller peaks in the Synchronization Scope,
% which do not turn the receiver on, reduce the synchronization threshold.

%%
close_system('zynqRadioWLAN80211BeaconRxAD9361AD9364SL',0);

%% References
% # IEEE Std 802.11-2007: _IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements, Part 11: Wireless LAN
% Medium Access Control (MAC) and Physical Layer (PHY) Specifications_, 
% IEEE, New York, NY, USA, 1999-2007.

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('zynqRadioWLAN80211BeaconRxAD9361AD9364SL_init.m')
% zynqRadioWLAN80211BeaconRxInitAD9361AD9364SL.m>
% * <matlab:edit('zynqRadioWLAN80211BeaconRxModelParamsAD9361AD9364SL.m')
% zynqRadioWLAN80211BeaconRxModelParamsAD9361AD9364SL.m>

displayEndOfDemoMessage(mfilename)

