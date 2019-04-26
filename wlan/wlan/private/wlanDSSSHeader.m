%wlanDSSSHeader DSSS PLCP Header
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   [Y, CRC] = wlanDSSSHeader(CFGDSSS) generates a DSSS modulated Physical
%   Layer Convergence Procedure (PLCP) Header time-domain waveform and PLCP
%   CRC field for the non-HT DSSS transmission format.
%
%   Y is the time-domain PLCP Header signal. It is a complex vector of
%   size Ns-by-1, where Ns represents the number of time-domain samples.
%
%   CRC is the PLCP CRC (CRC-16) field. It is a real vector of size
%   16-by-1.
%
%   CFGDSSS is the format configuration object of type <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>
%   which specifies the parameters for the Non-HT format. Only DSSS
%   modulation type is supported. DBPSK modulation is used for the long
%   preamble type and DQPSK modulation is used for the short preamble type;
%   CFGDSSS.DataRate does not affect the header modulation.
%
%   Example: 
%   %  Create the header waveform used for 802.11 DQPSK operation:
%   
%      cfgDSSS = wlanNonHTConfig('Modulation','DSSS');
%      cfgDSSS.DataRate = '2Mbps';
%
%      header = wlanDSSSHeader(cfgDSSS);
%
%   See also wlanNonHTConfig, wlanDSSSPreamble, wlanDSSSData.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

function [y,CRC] = wlanDSSSHeader(cfgDSSS)
    
    % Only applicable for DSSS modulation of the nonHT format configuration
    dsssValidateConfig(cfgDSSS,mfilename);
    % For DSSS DataRate='1Mbps', Preamble='Short' should be ignored and
    % Preamble='Long' used instead.
    if (strcmpi(cfgDSSS.Preamble,'Short') && strcmpi(cfgDSSS.DataRate,'1Mbps'))
        cfgDSSS.Preamble = 'Long';
    end
    
    % Create information structure
    cfgInfo = dsssInfo(cfgDSSS);

    % Create preamble bits; the preamble bits are scrambled and modulated
    % here to bring the scrambler and modulator to the correct state before
    % processing the header.
    % See dsssInfo(cfgDSSS) for details of SYNC and SFD fields
    preambleBits = [cfgInfo.Sync; cfgInfo.SFD];
    
    % Create CRC
    % See dsssInfo(cfgDSSS) for details of SIGNAL, SERVICE and LENGTH fields
    % Clause 17.2.3.7 PLCP CRC (CRC-16) field
    % (same as Clause 16.2.3.7 PLCP CRC field)
    % (same as Clause 17.2.3.14 Short CRC-16 field (shortCRC))
    CRC = dsssCRCGenerate([cfgInfo.Signal; cfgInfo.Service; cfgInfo.Length]);
    
    % Create header
    headerBits = [cfgInfo.Signal; cfgInfo.Service; cfgInfo.Length; CRC];
    
    % Scramble preamble + header
    % Clause 17.2.4 PLCP/High Rate PHY data scrambler
    % (same as Clause 16.2.4 PLCP/DSSS PHY data scrambler)
    scrambled = dsssScramble([preambleBits; headerBits], ...
                    cfgInfo.ScramblerInitialization);
    
    % Modulate preamble + header
    % 'P' is the position of the last preamble symbol in the symbol vector.
    P = length(preambleBits);
    if (strcmpi(cfgDSSS.Preamble,'Long'))
        % Clause 17.2.3.8 Long PLCP data modulation
        % (same as Clause 16.2.5 PLCP data modulation)
        % DBPSK modulation
        symbols = dsssPSKModulate(scrambled,'1Mbps');
    else
        % Clause 17.2.3.15 Short PLCP data modulation        
        % Repeat preamble to get DBPSK phase transitions through
        % the DQPSK modulator
        scrambled = [reshape(repmat(scrambled(1:P,1),1,2).',P*2,1); ...
                     scrambled(P+1:end,1)];
        % DQPSK modulation
        symbols = dsssPSKModulate(scrambled,'2Mbps');
    end
    
    % Extract modulated header
    symbols = symbols((P+1):end,1);
    
    % Spread modulated header
    % Clause 17.4.6.5 Spreading sequence and modulation for 1 Mb/s and 2 Mb/s
    % (same as Clause 16.4.6.4 Spreading sequence)
    % Barker spreading
    y = dsssBarkerSpread(symbols);
    
end

% [EOF]