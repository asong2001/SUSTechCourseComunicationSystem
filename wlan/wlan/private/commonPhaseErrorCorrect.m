function trackedSym = commonPhaseErrorCorrect(sym,cpe)
%commonPhaseErrorCorrect common phase error correction
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   TRACKEDSYM = commonPhaseErrorCorrect(SYM,CPE) returns phase corrected
%   OFDM symbols given a common phase error.
%
%   CPE is a 1-by-Nsym vector containing the common phase error per OFDM
%   symbol. Nsym is the number of OFDM symbols.
%
%   SYM is a complex Nsc-by-Nsym-by-Nr array containing the received OFDM
%   symbols at pilot subcarriers. Nsc is the number of subcarriers and Nr
%   is the number of receive antennas.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Phase correct
[Nsc,~,Nr] = size(sym);
x = exp(-1i*cpe); % Create constant
trackedSym = complex(zeros(size(sym)));
for r = 1:Nr
    for k = 1:Nsc
        trackedSym(k,:,r) = sym(k,:,r).*x;
    end
end
end