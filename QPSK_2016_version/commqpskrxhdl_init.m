function SimParams = commqpskrxhdl_init
% Set simulation parameters
% SimParams = sdruqpskrx_init

%   Copyright 2011-2012 The MathWorks, Inc.

% General simulation parameters
SimParams.M = 4; % M-PSK alphabet size
SimParams.Upsampling = 4; % Upsampling factor
SimParams.Downsampling = 2; % Downsampling factor
SimParams.Fs = 2e5; % Sample rate
SimParams.Ts = 1/SimParams.Fs; % Sample time
SimParams.FrameSize = 100; % Number of modulated symbols per frame

% Tx parameters
SimParams.BarkerLength = 13; % Number of Barker code symbols
SimParams.DataLength = (SimParams.FrameSize - SimParams.BarkerLength)*log2(SimParams.M); % Number of data payload bits per frame
SimParams.MsgLength = 105;
SimParams.RCGroupDelay = 5; % Group delay for the raised cosine receive filter

% Rx parameters
K = 1;
A = 1/sqrt(2);
% Look into model for details for details of PLL parameter choice.
SimParams.FineFreqPEDGain = 2*K*A^2+2*K*A^2; % K_p for Fine Frequency Compensation PLL, determined by 2KA^2 (for binary PAM), QPSK could be treated as two individual binary PAM
SimParams.FineFreqCompensateGain = 1; % K_0 for Fine Frequency Compensation PLL
SimParams.TimingRecTEDGain = 2.7*2*K*A^2+2.7*2*K*A^2; % K_p for Timing Recovery PLL, determined by 2KA^2*2.7 (for binary PAM), QPSK could be treated as two individual binary PAM, 2.7 is for raised cosine filter with roll-off factor 0.5
SimParams.TimingRecCompensateGain = -1; % K_0 for Timing Recovery PLL, fixed due to modulo-1 counter structure

load ('Data4commQPSKRxHDL', 'data_in');
assignin('base', 'cdfq', data_in);
load ('Data4commQPSKRxHDL', 'rcRxFilt1');
SimParams.rcRxFilt = rcRxFilt1;

% For CORDIC
SimParams.rad_WL = 16;
SimParams.rad_FL = 13;
SimParams.car_WL = 32;
SimParams.car_FL = 26;
SimParams.tb = [0.4636    0.2450    0.1244    0.0624    0.0312];

% The CORDIC output represents an angle, say alpha. Given M=3 and the fact 
% that CFC raises the input signal to the power of SimParams.M (4 in the 
% case of QPSK modulation), the actual frequency is estimated as
%
%    f_hat = alpha/[pi*SimParams.Ts*(3+1)]/SimParams.M;
%
% which, in one sample period, translates into a phase shift of
%
%    phase_hat = 2pi*f_hat*SimParams.Ts 
%              = alpha/(2*SimParams.M)
%
% This phase needs to multiply 2^17/(2pi) before feeding into NCO,
% therefore, the net constant should be (-1 is for compensation)
%
%    SimParams.CFC_Const = (-1)*[1/(2*SimParams.M)]*[2^17/(2pi)];
%                        = -(2^15)/(pi*SimParams.M)
% However, the (1/pi) term is already taken care of by Normalized Radians
% representation in Complex To Magnitude Angle HDL Optimized, so:

SimParams.CFC_Const = -(2^15)/(SimParams.M);

% Phase value needs to multiply 2^17/(2pi) before feeding into NCO,
% do not use HDl NCO block 
SimParams.FFC_NCOConst = -(2^16)/pi;

SimParams.WL_acc = 18;             % accumulator word length
SimParams.WL_lut = 16;             % w. length of quantized accumulator bits (input to LUT)
SimParams.WL_out=16;               % output word length
SimParams.FFC_Const = -2^(SimParams.WL_acc)/(2*pi);

SimParams.phMSB = SimParams.WL_acc-3;
SimParams.phLSB = SimParams.phMSB- SimParams.WL_lut+2+1;

invec = 0:2^(SimParams.WL_lut-2)-1;
SimParams.sin_tbl = sin(2*pi*invec./2^SimParams.WL_lut);
SimParams.cos_tbl = cos(2*pi*invec./2^SimParams.WL_lut);


