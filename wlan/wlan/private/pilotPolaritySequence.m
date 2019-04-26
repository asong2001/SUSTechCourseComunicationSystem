function seq = pilotPolaritySequence(varargin)
%pilotPolaritySequence Pilot polarity sequence
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   SEQ = pilotPolaritySequence(N) returns the pilot polarity sequence
%   specified in IEEE Std 802.11-2012 Eqn 18-25 for a given OFDM symbol
%   index N.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Pilot polarity sequence, p_n, IEEE Std 802.11-2012, eqn 18-25
pn = [1  1  1  1 -1 -1 -1 1 -1 -1 -1 -1  1  1 -1  1 -1 -1  1 1 -1  1  1 -1  1  1  1  1 ...
    1  1 -1  1  1  1 -1 1  1 -1 -1  1  1  1 -1  1 -1 -1 -1 1 -1  1 -1 -1  1 -1 -1  1 ...
    1  1  1  1 -1 -1  1 1 -1 -1  1 -1  1 -1  1  1 -1 -1 -1 1  1 -1 -1 -1 -1  1 -1 -1 ...
    1 -1  1  1  1  1 -1 1 -1  1 -1  1 -1 -1 -1 -1 -1  1 -1 1  1 -1  1 -1  1  1  1 -1 ...
    -1  1 -1 -1 -1  1  1 1 -1 -1 -1 -1 -1 -1 -1].';
NumPN = numel(pn);
if nargin==1
    idx =  varargin{1};   % Argument is OFDM symbol index (n)
else
    idx = 0:NumPN-1;  % Return whole sequence
end
% Cyclic extend sequence as required (modulus)
seq = pn(mod(idx,NumPN)+1);
end