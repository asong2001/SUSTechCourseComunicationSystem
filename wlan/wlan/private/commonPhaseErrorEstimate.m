function cpe = commonPhaseErrorEstimate(rxPilots,chanEstPilots,refPilots)
%commonPhaseErrorEstimate common phase error estimate
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   CPE = commonPhaseErrorEstimate(RXPILOTS,CHANESTPILOTS,REFPILOTS) returns the
%   common phase error per OFDM symbol, CPE. CPE is sized 1-by-Nsym, where
%   Nsym is the number of OFDM symbols.
%
%   RXPILOTS is a complex Nsp-by-Nsym-by-Nr array containing the received
%   OFDM symbols at pilot subcarriers. Nsp is the number of pilot
%   subcarriers and Nr is the number of receive antennas.
%
%   CHANESTPILOTS is a complex Nsp-by-Nsts-by-Nr array containing the
%   channel gains at pilot subcarriers. Nsts is the number of space-time
%   streams.
%
%   REFPILOTS is a complex Nsp-by-Nsym-by-Nsts array containing the
%   reference pilot values.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

[Np,Nsym,Nr] = size(rxPilots);
Nsts = size(chanEstPilots,2);

% Calculate an estimate of the received pilots using the channel estimate
temp = complex(zeros(Np,Nsym,Nsts,Nr));
for r = 1:Nr
    for s = 1:Nsts
        for k = 1:Np
            temp(k,:,s,r) = chanEstPilots(k,s,r).*refPilots(k,:,s);
        end
    end
end
% Sum over space-time streams and remove that dimension by permuting
estRxPilots = permute(sum(temp,3),[1 2 4 3]); 

% Phase correction based on Allert val Zelst and Tim C. W. Schenk,
% Implementation of a MIMO OFDM-Based Wireless LAN System, IEEE
% Transactions on Signal Processing, Vol. 52, No. 2, February 2004. The
% result is averaged over the number of receive antennas (summed over the
% 3rd dimension).
cpe = angle(sum(sum(rxPilots.*conj(estRxPilots),1),3));    
end