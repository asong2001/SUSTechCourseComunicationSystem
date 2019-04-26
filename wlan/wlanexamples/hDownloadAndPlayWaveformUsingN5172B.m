function hDownloadAndPlayWaveformUsingN5172B(instrumentAddress, ...
                        IQData, sampleRate, centerFrequency, outputPower)
% hDownloadAndPlayWaveformUsingN5172B Download and generate RF signal using
% N5172B.
%
%   hDownloadAndPlayWaveformUsingN5172B(instrumentAddress, ...
%    IQData, sampleRate, centerFrequency, outputPower) connects to an
%   Agilent N5172B Vector Signal Generator at the hostname/IP address
%   specified by instrumentAddress and downloads the complex baseband
%   signal specified by IQData and sampleRate to the instrument. It then
%   generates the RF signal at the specified centerFrequency and
%   outputPower. Note that sampleRate and centerFrequency need to be
%   specified in Hz and outputPower in dBm.
%
%   If the fcn call errors due to connection issues, retry after a 
%   clear all and instrreset.
%
%   Example: 
%    hDownloadAndPlayWaveformUsingN5172B('122.212.111.36', zeros(1,100), ...
%       1e6, 1.5e9, 0)

%   Copyright 2013-2014 The MathWorks Inc.

% Check for ICT presence
if isempty(ver('instrument')) 
    error('hDownloadAndPlayWaveformUsingN5172B:NoICT', ...
          ['Please install Instrument Control Toolbox for ', ...
           'this function to work.']);
end

%% Check input waveform and format it for download to instrument

% Check input waveform to ensure it is a vector with at least 60 points
if ~isvector(IQData) || numel(IQData)<60
    error('hDownloadAndPlayWaveformUsingN5172B:invalidInput','Invalid input data');
else
    IQsize = size(IQData);
    % User gave input as column vector. Reshape it to row vector.
    if ~isequal(IQsize(1),1)
        % warning('Wrong input detected. Automatically converting to row vector.');
        IQData = reshape(IQData,1,IQsize(1));
    end
end

% Separate out the real and imaginary data in the IQ Waveform
wave = [real(IQData);imag(IQData)];
wave = wave(:)'; % transpose the waveform

% Scale the waveform as necessary
tmp = max(abs([max(wave) min(wave)]));
if (tmp == 0)
    tmp = 1;
end

% ARB binary range is 2's Compliment -32768 to + 32767
% So scale the waveform to +/- 32767
scale = 2^15-1;
scale = scale/tmp;
wave = round(wave * scale);
modval = 2^16;
% Get data from double to uint16 as required by instrument
wave = uint16(mod(modval + wave, modval));

%% Download waveform to instrument

% Verify VISA installation and select VISA if available
foundVISA = instrhwinfo('visa');
if ~isempty(foundVISA.InstalledAdaptors)
    deviceObject = visa(foundVISA.InstalledAdaptors{1}, ...
                sprintf('TCPIP0::%s::inst0::INSTR',instrumentAddress));
    usingVISA = true;
else
    deviceObject = tcpip(instrumentAddress,5025);
    usingVISA = false;
end

% Set up the output buffer size
deviceObject.OutputBufferSize = 4*numel(wave) + 1024;
% Set the timeout
deviceObject.Timeout = 30;
% Set object to use BigEndian format
deviceObject.ByteOrder = 'bigEndian';
% filename for the data in the ARB
ArbFileName = 'MATLABWfm';

% Open connection to the instrument
fopen(deviceObject);

% Clear hardware buffers on the instrument
if usingVISA
    % This is only possible with VISA objects
    clrdevice(deviceObject);
end

% Reset instrument and clear instrument error queue
fprintf(deviceObject,'*RST;*CLS');

% Some settings commands to make sure we don't damage the instrument
fprintf(deviceObject,':SOURce:RADio:ARB:STATe OFF');
fprintf(deviceObject,':OUTPut:MODulation:STATe OFF');
fprintf(deviceObject,':OUTPut:STATe OFF');

% Write the data to the instrument
fprintf('Starting Download of %d IQ samples...\n',size(wave,2)/2);
binblockwrite(deviceObject,wave,'uint16',[':MEM:DATA "NVWFM:' ArbFileName '",']);
fprintf(deviceObject,'');
% Wait till operation completes
localWaitTillComplete;

% Clear volatile memory waveforms
fprintf(deviceObject, ':MMEMory:DELete:WFM');
% Copy the waveform to volatile memory
fprintf(deviceObject,[':MEMory:COPY:NAME "NVWFM:' ArbFileName '","WFM1:NVWFM"']);
% Wait till operation completes
localWaitTillComplete;
% Display any instrument errors
localDisplayInstrumentErrors;

fprintf('\t...done!\n');

% Set up markers if we need to use this for synchronization using the 
% Event1 hardware output on the signal generator

% Clear all markers from the file
fprintf(deviceObject,':SOURce:RADio:ARB:MARKer:CLEar:ALL "NVWFM",1');
% Set marker 1 (first input after filename), starting at the first point
% (second input), ending at point 1 (third input) and skipping 0.
% Refer page 295 of 
% <http://cp.literature.agilent.com/litweb/pdf/N5180-90004.pdf Programmer's manual> 
% for more info
fprintf(deviceObject,':SOURce:RADio:ARB:MARKer:SET "NVWFM",1,1,1,0');
% Play back the selected waveform
fprintf(deviceObject, ':SOURce:RADio:ARB:WAVeform "WFM1:NVWFM"');

%% Generate the signal
% Display message
fprintf('Generating RF Signal...\n');
% Set the sample rate (Hz) for the signal.
% You can get this info for the standard signals by looking at the data 
% in the 'waveforms' variable
fprintf(deviceObject,[':SOURce:RADio:ARB:SCLock:RATE ' num2str(sampleRate) 'Hz']);
% set center frequency (Hz)
fprintf(deviceObject, ['SOURce:FREQuency ' num2str(centerFrequency) 'Hz']);
% set output power (dBm)
fprintf(deviceObject, ['POWer ' num2str(outputPower)]);
% set runtime scaling to 75% of DAC range so that interpolation by the
% instrument baseband generator doesn't cause errors
% Refer page 159 of
% <http://cp.literature.agilent.com/litweb/pdf/E4400-90503.pdf User guide> 
% for more info
fprintf(deviceObject, 'RADio:ARB:RSCaling 75');
% Turn on output protection 
fprintf(deviceObject,':OUTPut:PROTection ON');
% ARB Radio on
fprintf(deviceObject, ':SOURce:RADio:ARB:STATe ON');
% modulator on
fprintf(deviceObject, ':OUTPut:MODulation:STATe ON');
% RF output on
fprintf(deviceObject, ':OUTPut:STATe ON');
% Display any instrument errors
localDisplayInstrumentErrors;
fprintf('\t...done!\n');

%% Clean up

% Close instrument connection and delete object
fclose(deviceObject); 
delete(deviceObject); 

    function localDisplayInstrumentErrors()
    % Display any instrument errors
        instrumentError = query(deviceObject,'SYSTem:ERRor?');
        while isempty(strfind(lower(instrumentError),'no error'))
            fprintf('\tInstrument Error: %s',instrumentError);
            instrumentError = query(deviceObject,'SYSTem:ERRor:NEXT?');
        end
    end

    function localWaitTillComplete()
    % Wait until instrument operation is complete
        operationComplete = str2double(query(deviceObject,'*OPC?'));
        while ~operationComplete
            operationComplete = str2double(query(deviceObject,'*OPC?'));
        end
    end
end