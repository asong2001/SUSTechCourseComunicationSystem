function out = wlanCyclicShift(in, cSh, Nfft, mode, varargin)
%wlanCyclicShift Cyclic shift delay insertion in frequency domain
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   out = wlanCyclicShift(IN, CSH, NFFT, MODE) performs cyclic shift delay
%   insertion to all subcarriers in the frequency domain signal IN, where
%   size(IN,2)=length(CSH). The other function parameters are
%       CSH is the shift value in samples (from getCyclicShiftVal)
%       MODE is 'Tx' or 'Rx'
%       NFFT is the length of the FFT
%
%   out = wlanCyclicShift(IN, CSH, NFFT, MODE, INDICES) performs cyclic
%   shift delay insertion in the subcarrier locations specified in INDICES.
%
%   See also getCyclicShiftVal.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

if nargin>4
    ind = varargin{1};
else
    ind = 1:Nfft;
end

% Invert shift for receive processing
if strcmp(mode,'Rx')
    cSh = -1 * cSh;
end

% Apply cyclic shift
k = (1:Nfft) - Nfft/2 -1;
phaseShift = exp(-1i*2*pi*cSh*k/Nfft).';
out = bsxfun(@times,in,phaseShift(ind,:));

end

% [EOF]