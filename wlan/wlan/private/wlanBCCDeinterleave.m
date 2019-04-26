function y = wlanBCCDeinterleave(x, Format, Ncbps, Nbpscs, varargin)
%wlanBCCDeinterleave BCC Deinterleaver
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   y = wlanBCCDeinterleave(X, FORMAT, NCBPS, NBPSCS, CBW, NSS) outputs
%   the deinterleaved input (X) using the following parameters
%     FORMAT: either of VHT, HT_MF or NON_HT
%     NCBPS : Input block length per stream (Ncbps or Ncbpssi)
%     NBPSCS: Modulation order
%     CBW   : Channel bandwidth per segment (one of 20, 40, 80)
%     NSS   : Number of spatial streams
%
%   For the NON_HT format, CBW and NSS are not required parameters.
%
%   wlanBCCDeinterleave(randi([0 1], 52, 1), 'VHT', 52, 1, 20, 1,)
%
%   See also wlanBCCInterleave.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(4,6);
assert(isreal(x));

s = max(Nbpscs/2, 1);
if strcmp(Format, 'NON_HT')
    pNcol  = 16;                   % fixed value
    pNrow  = Ncbps/pNcol;
    
    % first stage deinterleaving - Eq. 22-82
    piElem = ( s*floor( (0:Ncbps-1)/s ) + ...
        mod( ((0:Ncbps-1) + floor( pNcol*(0:Ncbps-1)/Ncbps ) ), s ) ...
        + 1).';   % 1-based indexing for interleaving

else % HT and VHT
    CBW = varargin{1};
    Nss = varargin{2}; % should also be equal to size(inp,2)

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
        
    % Second stage deinterleaver - Eq. 22-82
    piElem = ( s*floor( (0:Ncbps-1)/s ) + ...
        mod( ((0:Ncbps-1) + floor( pNcol*(0:Ncbps-1)/Ncbps ) ), s ) ...
        + 1).';   % 1-based indexing for interleaving
    
    % First, if needed, permutation tables for Nss > 1
    pRMat = zeros(Ncbps, Nss);
    pRMat(:,1) = (1:Ncbps).';
    if Nss >= 2 && Nss <=4
        % Eq. 22-80, pg 282
        for iss = 2:Nss
            pRMat(:, iss) = mod((0:Ncbps-1).' + ( mod(2*(iss-1),3) + ...
                3*floor((iss-1)/3)) * pNrot * Nbpscs, Ncbps) + 1;
        end
    else
        jTab = [0 5 2 7 3 6 1 4];   % Table 22-18
        % Eq 22-81, pg 282
        for iss = 2:Nss
            pRMat(:, iss) = mod((0:Ncbps-1).' + jTab(iss) * ...
                pNrot * Nbpscs, Ncbps) + 1;
        end
        
    end
end

% Step
y = zeros(size(x)); inp1 = zeros(size(x,1),1);
if strcmp(Format, 'NON_HT')
    inp1(piElem, 1) = x;                    % second
    tmp = reshape(inp1, pNrow, pNcol).';    % first
    y = tmp(:);

else % HT or VHT
    inp2 = inp1;
    for nssIdx = 1:Nss
        inp1(pRMat(:,nssIdx), 1) = x(:, nssIdx);    % third
        inp2(piElem, 1) = inp1;                     % second        
        tmp = reshape(inp2, pNrow, pNcol).';        % first
        y(:, nssIdx) = tmp(:);
    end
    
end

end

% [EOF]
