%% 802.11a(TM) Transmission and Reception Using Analog Devices AD9361/AD9364
% This example shows how to use the Xilinx(R) Zynq-based radio support
% package with MATLAB(R) and WLAN System Toolbox(TM) to generate a
% simultaneous transmission and reception on a single SDR platform.
%
%% Example Summary
% WLAN System Toolbox can be used to generate standard-compliant baseband
% IQ waveforms. These baseband waveforms can be modulated for RF
% transmission using SDR radio hardware such as Xilinx Zynq-based radio.
%
% In this example, an image file is imported and packed into multiple WLAN
% packets of a baseband waveform that is generated using WLAN System
% Toolbox. A single antenna is used to generate an IEEE(R) 802.11a(TM)
% waveform. The RF WLAN waveform is then created, the baseband waveform is
% transferred to the hardware memory on the Zynq radio and transmitted over
% the air.
%
% The RF card used in this example is capable of simultaneous transmission
% and reception. Therefore, the transmitted signal is captured using the
% same Zynq radio hardware platform. The diagram below shows the setup
% used.
%
% <<SDRWLANSISO80211aTransceiverZynq_published.png>>
%
% The receiver captures a number of WLAN packets and performs
% synchronization, channel estimation and equalization to retrieve packet
% parameters. The data field is then extracted and the transmitted payload
% is recovered using the retrieved packet parameters. After decoding the
% received waveform, the transmitted image is recovered.
%
%% SDR Support Packages
% This example requires the Xilinx Zynq-based radio support package. This
% can be installed using the <matlab:supportPackageInstaller support
% package installer>.
%
% More information about other supported SDR platforms can be found <http://www.mathworks.com/hardware-support/index.html?q=%20product:%22Communications+System+Toolbox%22
% here>.
%
%% Full Example 
% The full example description and source code can be found in the 
% <http://www.mathworks.com/help/supportpkg/xilinxzynqbasedradio/examples.html?refresh=true list of examples using Xilinx Zynq-Based Radio>
% under the name "802.11a Transmission and Reception Using Analog Devices
% AD9361/AD9364".

% Copyright 2015 The MathWorks, Inc.

displayEndOfDemoMessage(mfilename)