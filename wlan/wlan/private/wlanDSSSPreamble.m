function y = wlanDSSSPreamble(cfgDSSS)
%wlanDSSSPreamble DSSS PLCP Preamble
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanDSSSPreamble(CFGDSSS) generates a DSSS modulated PLCP Preamble
%   time-domain waveform for the non-HT DSSS transmission format.
%
%   Y is the time-domain PLCP Preamble signal. It is a complex vector of
%   size Ns-by-1, where Ns represents the number of time-domain samples.
%
%   CFGDSSS is the format configuration object of type <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>
%   which specifies the parameters for the Non-HT format. Only DSSS
%   modulation type is supported. DBPSK modulation is used for the
%   preamble; CFGDSSS.DataRate does not affect the preamble modulation. 
%
%   Example: 
%   %  Create the short preamble waveform used for 802.11b HR/DSSS/short 
%   %  operation.
%   
%      cfgDSSS = wlanNonHTConfig('Modulation','DSSS');
%      cfgDSSS.DataRate = '11Mbps';
%      cfgDSSS.Preamble = 'Short';
%
%      preamble = wlanDSSSPreamble(cfgDSSS);
%
%   See also wlanNonHTConfig, wlanDSSSHeader, wlanDSSSData.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

    % Only applicable for DSSS non-HT modulation configuration
    dsssValidateConfig(cfgDSSS,mfilename);
    % For DSSS DataRate='1Mbps', Preamble='Short' should be ignored and
    % Preamble='Long' used instead.
    if (strcmpi(cfgDSSS.Preamble,'Short') && strcmpi(cfgDSSS.DataRate,'1Mbps'))
        cfgDSSS.Preamble = 'Long';
    end
    
    % Create information structure
    cfgInfo = dsssInfo(cfgDSSS);
    
    % Concatenate SYNC and SFD fields
    % See dsssInfo(cfgDSSS) for details
    preambleBits = [cfgInfo.Sync; cfgInfo.SFD];
    
    % Scramble preamble
    % Clause 17.2.4 PLCP/High Rate PHY data scrambler
    % (same as Clause 16.2.4 PLCP/DSSS PHY data scrambler)
    scrambled = dsssScramble(preambleBits,cfgInfo.ScramblerInitialization);
    
    % Modulate preamble
    % Clause 17.2.3.8 Long PLCP data modulation
    % (same as Clause 17.2.3.15 Short PLCP data modulation)
    % (same as Clause 16.2.5 PLCP data modulation)
    % DBPSK modulation
    symbols = dsssPSKModulate(scrambled,'1Mbps');
    
    % Spread preamble
    % Clause 17.4.6.5 Spreading sequence and modulation for 1 Mb/s (and 2Mb/s)
    % (same as Clause 16.4.6.4 Spreading sequence)
    % Barker spreading
    y = dsssBarkerSpread(symbols);
    
end

% [EOF]