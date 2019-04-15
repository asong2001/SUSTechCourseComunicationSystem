clc;
clear;

% call commQPSKTXRX simulation
SNR = 0:1:20;

for k = 1:length(SNR)
    clear all;
    EbNo = SNR(k);
    ber(k) = commQPSKTransmitterReceiver(EbNo);
end