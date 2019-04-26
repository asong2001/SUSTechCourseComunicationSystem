%% 802.11(TM) OFDM Beacon Receiver with USRP(R) Hardware
% This example shows how to use the Universal Software Radio Peripheral
% (USRP(R)) device using SDRu (Software Defined Radio USRP) System
% objects(TM) to implement a WLAN receiver. The receiver is able to recover
% 802.11(TM) OFDM Non-HT beacon frames transmitted over the air from
% commercial 802.11 hardware.
%
%% Example Summary
% WLAN System Toolbox(TM) provides functions and tools to decode 802.11
% waveforms. This example shows how to receive signals from commercial WLAN
% transmitters in MATLAB using USRP. A receiver design is demonstrated
% including synchronization, transmission configuration recovery, and
% payload decoding for Non-HT packets.
%
% <<usrpBeaconSDRDiagram.png>>
%
% In this example, OFDM beacon packets corrupted by the transmission over
% the air are captured and processed to recover the payload contents. To
% recover beacon packets the receiver performs packet detection, symbol
% timing and frequency offset correction, and L-SIG and payload decoding.
% The resulting payload bits are then evaluated to determine whether the
% payload is a beacon frame and the contents are displayed as appropriate.
%
% To view an example of a WLAN front end processing which does not require
% SDR hardware see <NonHTBeaconReceiverExample.html 802.11(TM) OFDM Beacon
% Receiver with Live Data>.

%% SDR Support Packages
% This example requires the USRP-based radio support package. This can be
% installed using the <matlab:supportPackageInstaller support package
% installer>.
%
% More information about other supported SDR platforms can be found <http://www.mathworks.com/hardware-support/index.html?q=%20product:%22Communications+System+Toolbox%22
% here>.
%
%% Full Example 
% The full example description and source code can be found in the 
% <http://www.mathworks.com/help/supportpkg/usrpradio/examples.html?refresh=true list of examples using USRP>
% under the name "IEEE 802.11(TM) WLAN - OFDM Beacon Receiver with USRP(R) Hardware".

% Copyright 2015 The MathWorks, Inc.

displayEndOfDemoMessage(mfilename)