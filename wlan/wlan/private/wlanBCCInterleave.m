function y = wlanBCCInterleave(x, Format, Ncbps, Nbpscs, varargin)
%wlanBCCInterleave BCC Interleaver.
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   y = wlanBCCInterleave(X, FORMAT, NCBPS, NBPSCS, CBW, NSS) outputs
%   the interleaved binary convolutionally encoded data input (X) using
%   the following parameters
%     FORMAT: either of VHT, HT_MF or NON_HT
%     NCBPS : Input block length per stream (Ncbps or Ncbpssi)
%     NBPSCS: Modulation order
%     CBW   : Channel bandwidth per segment (one of 20, 40, 80)
%     NSS   : Number of spatial streams
%
%   y = wlanBCCInterleave(X, FORMAT, NCBPS, NBPSCS) is used for the NON_HT
%   format, CBW and NSS are not required parameters.
%
%   Example:
%   % Perform VHT interleaving.
%   inBits = randi([0 1], 52, 1);
%   Ncbps = length(inBits);
%   Nbpscs = 1;
%   CBW = 20;
%   Nss = 1;
%   out = wlanBCCInterleave(inBits, 'VHT', Ncbps, Nbpscs, CBW, Nss);
%
%   See also wlanBCCDeinterleave.

%   Copyright 2015 The MathWorks, Inc.


%#codegen

% Validate inputs
narginchk(4,6);

s = max(Nbpscs/2, 1);   % Eq. 22-68, IEEE Std 802.11ac-2013
if strcmp(Format, 'NON_HT')
    pNcol  = 16;  % fixed value, Eq. 18-18
    pNrow  = Ncbps/pNcol;
    
    % For second stage interleaver - Eq. 22-77
    piElem = ( s*floor( (0:Ncbps-1)/s ) + ...
        mod( ((0:Ncbps-1) + Ncbps - ...
        floor( pNcol*(0:Ncbps-1)/Ncbps ) ), s ) ...
        + 1).';   % 1-based indexing for interleaving
else % HT and VHT
    CBW = varargin{1};
    Nss = varargin{2}; % should also be equal to size(inp,2)

    % Rows and columns in interleaver, Table 22-17, IEEE Std 802.11ac-2013
    switch CBW
        case {20}
            pNcol = 13;
            pNrow = 4*Nbpscs;
            if Nss <= 4
                pNrot = 11;
            else
                pNrot = 6;
            end
        case {40}
            pNcol = 18;
            pNrow = 6*Nbpscs;
            if Nss <= 4
                pNrot = 29;
            else
                pNrot = 13;
            end
        otherwise % {80}
            pNcol = 26;
            pNrow = 9*Nbpscs;
            if Nss <= 4
                pNrot = 58;
            else
                pNrot = 28;
            end
    end
    
    % For second stage interleaver  - Eq. 22-77
    piElem = ( s*floor( (0:Ncbps-1)/s ) + ...
        mod( ((0:Ncbps-1) + Ncbps - ...
        floor( pNcol*(0:Ncbps-1)/Ncbps ) ), s ) ...
        + 1).';   % 1-based indexing for interleaving
    
    % Third, if needed, permutation tables for Nss > 1
    pRMat = zeros(Ncbps, Nss);
    pRMat(:,1) = (1:Ncbps).';
    if Nss >= 2 && Nss <=4
        % Eq. 22-78, pg 281
        for iss = 2:Nss
            pRMat(:, iss) = mod((0:Ncbps-1).' - ( mod(2*(iss-1),3) + ...
                3*floor((iss-1)/3)) * pNrot * Nbpscs, Ncbps) + 1;
            % 1-based indexing for interleaving
        end
    else % Applies for Nss > 4
        jTab = [0 5 2 7 3 6 1 4];   % Table 22-18
        % Eq. 22-79
        for iss = 2:Nss
            pRMat(:, iss) = mod((0:Ncbps-1).' - jTab(iss) * ...
                pNrot * Nbpscs, Ncbps) + 1; % 1-based indexing for interleaving
        end
        
    end
end

% Interleave input data
y = zeros(size(x));
if strcmp(Format, 'NON_HT')
    tmp = reshape(x, pNcol, pNrow).';
    y(piElem, :) = tmp(:);
    
else % HT or VHT
    
    inp2 = zeros(size(x,1),1);
    for nssIdx = 1:Nss
        tmp = reshape(x(:, nssIdx), pNcol, pNrow).';
        inp2(piElem, :) = tmp(:);
        y(pRMat(:,nssIdx), nssIdx) = inp2;
    end
end

end

% [EOF]
