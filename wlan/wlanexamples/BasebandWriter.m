classdef BasebandWriter < matlab.System
    %BasebandWriter Write baseband data to a compressed file.
    %   S = BasebandWriter returns a System object S that writes
    %   complex-valued data to a compressed data file.  One frame of data
    %   is written at a time.
    %
    %   S = BasebandWriter('PropertyName',PropertyValue, ...) returns a
    %   System object with named properties set to specified values.
    %
    %   S = BasebandWriter(FNAME,FS,FC,MaxS,'PropertyName',PropertyValue,
    %   ...) returns a BasebandWriter object where 'Filename' is set to
    %   FNAME, 'SampleRate' is set to FS, 'CenterFrequency' is set to FC,
    %   and 'MaxSamples' is set to MaxS.  If a file extension is omitted
    %   from FNAME, '.bb' is assumed.  These values can also be set using
    %   property name and value pairs.
    %
    %   A baseband data file is a specific type of MAT file written by
    %   BasebandWriter. Baseband data is a matrix of complex-valued (I/Q)
    %   samples of a signal that have been down-converted from a center
    %   frequency of FC Hz to 0 Hz ("baseband"), and sampled at a rate of
    %   FS Hz. The matrix contains consecutive time samples down each
    %   column, and one column per data channel. FC and FS are recorded
    %   when a baseband data file is created.  Additional meta-data may be
    %   stored in the UserData property and is written to the baseband
    %   file.
    %
    %   Samples are written to the BasebandWriter object when the step
    %   method is called:
    %       step(S,BB);
    %   BB is a SamplesPerFrame-by-NumChannels matrix of complex values.
    %
    %   When the object S is cleared from memory or goes out of scope, any
    %   data remaining in S is written to the baseband file.  Call close(S)
    %   or release(S) to force this action without clearing S.  Calling
    %   close(S) allows subsequent writes to occur to the object, while
    %   release(S) deletes the object and no further writes may occur.
    %
    %   See also BasebandReader, convertFile, convertData.

    % Copyright 2015 The MathWorks, Inc.
    
    properties (Nontunable,Dependent,Hidden)
        % Number of SAMPLES PER CHANNEL to write to data file at a time.
        % Larger values increase performance at the expense of increased
        % memory usage.
        CacheSize
    end
    
    properties (Nontunable,Hidden)
        % Set to 'manual' to set cache size to value of CacheSize
        % property, or set to 'auto' to automatically set cache size based
        % on file length.
        %
        % If value of CacheSize property is modified, CacheSizeMode
        % automatically changes to 'manual'.
        CacheSizeMode = 'auto'
    end
    
    properties (Nontunable)
        % Name of baseband data file.  If a file extension is omitted,
        % '.mat' is assumed.
        Filename = 'data'
        
        % Center frequency of signal prior to downconversion to baseband,
        % in Hz.
        CenterFrequency = 1e9
        
        % Sample rate of baseband samples, in Hz.
        SampleRate = 200e3
        
        % Maximum number of most recent samples to retain in file.
        % Set to inf to record all samples.
        MaxSamples = inf
    end
    
    properties (SetAccess=private,Dependent)
        % Index of last sample written to file.
        SamplesWritten
    end
    
    properties (Nontunable)
        % User-specified data stored with the object.
        UserData
    end
    
    properties (Access=private)
        % Struct of internal properties, for efficiency
        %   .pMatObj
        %   .pFilename
        %   .pStepFcn
        %   .pLastSampleWritten
        %   .pMaxSampleWritten
        %
        %   .pCacheSampleRange - maps absolute time index to cache samples, [firstIdx lastIdx]
        %   .pCacheSamples     - cache samples / data values
        %   .pSamplesInCache   - # samples that can fit in cache
        %   .pNeedFlush        - T/F
        pInternal
        
        pCacheSize = 0
    end
    
    methods
        function sys = BasebandWriter(varargin)
            setProperties(sys,nargin,varargin{:}, ...
                'Filename','SampleRate','CenterFrequency', ...
                'MaxSamples');
        end
        
        function set.Filename(sys,val)
            if ~ischar(val)
                error('Filename must be a string.');
            end
            sys.Filename = val;
            
            % Reset internal properties:
            sys.pInternal = []; %#ok<MCSUP>
        end
        
        function set.CenterFrequency(sys,val)
            % Center frequency can be positive, zero or negative frequency,
            % since the samples are generally complex.
            validateattributes(val,{'numeric'}, ...
                {'scalar','real','finite'}, ...
                mfilename('class')','CenterFrequency');
            sys.CenterFrequency = val;
        end
        
        function set.MaxSamples(sys,val)
            validateattributes(val,{'numeric'}, ...
                {'scalar','real','positive'}, ...
                mfilename('class')','MaxSamples');
            sys.MaxSamples = val;
        end
        
        function set.SampleRate(sys,val)
            validateattributes(val,{'numeric'}, ...
                {'scalar','real','positive','finite'}, ...
                mfilename('class')','SampleRate');
            sys.SampleRate = val;
        end
        
        function y = get.CacheSize(sys)
            y = sys.pCacheSize;
        end
        
        function set.CacheSize(sys,val)
            validateattributes(val,{'numeric'}, ...
                {'scalar','real','positive','integer','finite'}, ...
                mfilename('class'),'CacheSize');
            sys.pCacheSize = val;
            
            % Side-effect:
            %
            % This will reset pInternal as well, which is desired.
            sys.CacheSizeMode = 'manual';
        end
        
        function set.CacheSizeMode(sys,val)
            sys.CacheSizeMode = validatestring(val, ...
                {'auto','manual'}, mfilename('class'),'CacheSizeMode');
            
            sys.pInternal = []; %#ok<MCSUP>
        end
        
        function y = get.SamplesWritten(sys)
            p = sys.pInternal;
            if isempty(p)
                y = [];
            else
                y = p.pLastSampleWritten;
            end
        end
        
        function close(sys)
            % Finalize MAT file
            p = sys.pInternal;
            if ~isempty(p) && p.pNeedFlush
                flushCache(sys); % must come before linearization
                linearizeCircularData(sys);
                
                % Clear flag so flush/linearize doesn't reoccur
                p.pNeedFlush = false;
                sys.pInternal = p;
            end
        end
    end
    
    methods (Static)
        function convertFile(varargin)
            %convertFile Convert legacy file into a baseband file.
            %  BasebandWriter.convertFile('LegacyFile','VarName', ...
            %  'NewFile',FS,FC) reads data from variable VarName stored in
            %  legacy MAT file named LegacyFile, and writes it into a
            %  baseband data file named NewFile with sample rate FS and
            %  center frequency FC. FS and FC take default values if
            %  omitted.
            %
            %  BasebandWriter.convertFile('LegacyFile','VarName', ...
            %  'NewFile','FS','FC') reads sample rate from variable 'FS'
            %  and center frequency from variable 'FC' in LegacyFile.
            
            % LegacyFile specified
            % convertData(LegacyFile,VarName,NewFile)
            % convertData(LegacyFile,VarName,NewFile,FS)
            % convertData(LegacyFile,VarName,NewFile,FS,FC)
            narginchk(3,5);
            
            % Option 1: if legacy file is saved using -v7.3 or later, we
            % can stream the data out of the file
            %sys = BasebandWriter(varargin{3:end});
            %readMatFile(sys,varargin{1:2});
            
            % Option 2: if legacy file is saved using -v7 or earlier, we
            % need to read it all into memory at once
            s = load(varargin{1},'-mat',varargin{2});
            data = s.(varargin{2}); % get named variable from struct
            % 'FS' and 'FC' may be specified instead of values.  Resolve
            % them from the named variables read from the file:
            if nargin > 3 && ischar(varargin{4})
                varargin{4} = s.(varargin{4});
            end
            if nargin > 4 && ischar(varargin{5})
                varargin{5} = s.(varargin{5});
            end
            BasebandWriter.convertData(data,varargin{3:end});
        end
        
        function convertData(varargin)
            %convertData Convert workspace data into a baseband file.
            %  BasebandWriter.convertData(X,NewFile,Fs,Fc) stores
            %  matrix X in the MATLAB workspace, with one column per
            %  channel and time samples down each column, into a baseband
            %  data file named NewFile with sample rate Fs and center
            %  frequency Fc. Fs and Fc take default values if omitted.
            
            % Data specified
            % convertData(X,NewFile)
            % convertData(X,NewFile,Fs)
            % convertData(X,NewFile,Fs,Fc)
            narginchk(2,4);
            sys = BasebandWriter(varargin{2:end});
            step(sys,varargin{1});
            
            % sys will go out of scope so cache is guaranteed to flush.
            % We call release/close explicitly nevertheless:
            release(sys);
        end
    end
    
    methods (Access=protected)
        function setupImpl(sys,X)
            % Create pInternal structure, providing fast property updates.
            % Reset cache explicitly.
            
            % Choose MAT-file API:
            %   false: public API, slower file-writes
            %    true: internal-only API, faster file-writes
            p.pFastWrite = false;
            p.pNeedFlush = false;
            
            p = checkInputData(sys,X,p);
            p = checkFileNameAndDeleteFile(sys,p);
            p = initMatFileObject(sys,p);
            p = chooseStepFcn(sys,p);
            p = initCacheSize(sys,p);
            p = initFastWrite(sys,p);
            
            % Use struct of private properties for access speed
            sys.pInternal = p;
            
            resetCacheAndFile(sys);
        end
        
        function resetImpl(sys)
            resetCacheAndFile(sys);
        end
        
        function stepImpl(sys,X)
            % Write samples to baseband data file.
            feval(sys.pInternal.pStepFcn,sys,X);
        end
        
        function releaseImpl(sys)
            close(sys);
        end
        
        function ret = infoImpl(sys)
            
            p = sys.pInternal;
            if isempty(p)
                s.Filename = '';
                s.ComplexData = [];
                s.DataClass = '';
                s.NumChannels = [];
                s.SamplesPerFrame = [];
                s.BytesInCache = [];
                s.SamplesInCache = [];
                s.FramesInCache = [];
            else
                s.Filename = p.pFilename;
                scalarVal = p.pLikeValue;
                s.ComplexData = ~isreal(scalarVal);
                s.DataClass = class(scalarVal);
                s.NumChannels = p.pNumChannels;
                s.SamplesPerFrame = p.pSamplesPerFrame;
                s.BytesInCache = p.pSamplesInCache ...
                    * bytesPerSample(scalarVal);
                s.SamplesInCache = p.pSamplesInCache;
                s.FramesInCache = p.pSamplesInCache / p.pSamplesPerFrame;
            end
            
            if nargout > 0
                ret = s;
            else
                disp(s);
            end
        end
    end
    
    methods (Access=private)
        function readMatFile(obj,F,V)
            %readMatFile Add data from a variable in a MAT file.
            %   readMatFile(SYS,F,V) stores data from variable V in MAT
            %   file F to the baseband data file object SYS. Unless file F
            %   was saved using MAT file version '-v7.3' or later, all data
            %   will be read into memory before it is stored into SYS.
            %
            %   Variable V must be numeric matrix of size MxN, with M time
            %   samples per channel and N channels.  After conversion, use
            %   BasebandDataReader to efficiently stream data from the new
            %   baseband data file.
            %
            %   Example: Read data from variable 'y' in handel.mat
            %            and add it to a new baseband data file.
            %   Fc = 0;     % audio has 0 Hz center frequency
            %   Fs = 8192;  % audio sample rate
            %   SYS = BasebandWriter('handelStream',Fc,Fs)
            %   readMatFile(SYS,'handel.mat','y');
            
            validateattributes(F,{'char'}, ...
                {'vector','row'}, ...
                mfilename('class'),'FileName',2);
            
            validateattributes(V,{'char'}, ...
                {'vector','row'}, ...
                mfilename('class'),'VariableName',3);
            
            oldMat = matfile(F);
            Vr = size(oldMat,V,1); %#ok<GTARG> % # rows in variable
            
            warnState = warning('off','MATLAB:MatFile:OlderFormat');
            
            % Get cache size.
            % If there were no prior calls to step(), it will be zero.
            % In that case, call initCacheSize early.
            %
            Nc = obj.pCacheSize;  % # samples per frame for BDW
            if Nc==0
                scalarValue = oldMat.(V)(1,1);
                obj.pInternal = initCacheSize(obj,obj.pInternal,scalarValue);
                Nc = obj.pCacheSize;
            end
            
            % Try to stream data from legacy MAT file.  If it's not
            % '-v7.3', user will get a warning -- but it's going to do the
            % job.
            i1 = 1;
            while i1 <= Vr
                i2 = min(Vr,i1+Nc-1); % stop at end of variable
                step(obj,oldMat.(V)(i1:i2,:));
                i1 = i2+1;
            end
            % don't close object, leave that to caller
            
            warning(warnState);
        end

        function p = initFastWrite(sys,p) %#ok<INUSL>
            % For more efficient MAT file reads, we use:
            %    y = matlab.internal.language.partialLoad(src,var,'-mat')
            % where
            %   src: fully-qualified MAT file name
            %   var: struct with the fields:
            %      .Name = 'data'
            %      .Subsets.Type = '()'
            %              .Index = {[sampleFirst sampleStep sampleLast], ...
            %                        [chanFirst chanStep chanLast]}
            % sampleStep can be omitted.
            % Dummy values are recorded for sampleFirst and sampleLast.
            %
            % This requires removal of the private-property
            % attribute from the following classes:
            %   matlab.internal.language.Subset: .Index
            %   matlab.internal.language.VariableSubset: .Subsets
            if p.pFastWrite
                subset = matlab.internal.language.Subset('()', ...
                    {[1 1],[1 p.pNumChannels]});
                p.pPartialSaveStruct = ...
                    matlab.internal.language.VariableSubset('data',subset);
            end
        end

        function p = initCacheSize(sys,p,scalarValue)
            % Init cache size
            %
            % Requires p.pSamplesPerFrame
            
            if nargin < 3
                scalarValue = p.pLikeValue;
                isFullDataSet = false;
            else
                isFullDataSet = true;
            end
            
            % Determine max # samples that fit in 256kB:
            bigS = floor(256*1024 / bytesPerSample(scalarValue));
            maxS = sys.MaxSamples;
            if isFullDataSet
                % No explicit frame size given, since it's a full dataset
                % (e.g., importing a full legacy MAT file).
                %
                % This works even if isinf(maxS)
                c = min(maxS,bigS);
            elseif strcmpi(sys.CacheSizeMode,'auto')
                % This works even if isinf(maxS)
                c = max(p.pSamplesPerFrame,min(maxS,bigS));
            else
                % Manual cache size specified.
                %
                % Use the greater of manual spec and spf:
                c = max(p.pSamplesPerFrame,sys.pCacheSize);
            end
            sys.pCacheSize = c;
            p.pSamplesInCache = c; % actual # cache samples per channel
        end
        
        function p = chooseStepFcn(sys,p)
            % Choose the step function.
            % Needs p.pSamplesPerFrame
            
            maxS = sys.MaxSamples;
            if isinf(maxS)
                % Unlimited file size.
                % Always append data to end.
                p.pStepFcn = @writeData;
            else
                if p.pSamplesPerFrame <= maxS
                    p.pStepFcn = @stepMaxShortData;
                else
                    p.pStepFcn = @stepMaxLongData;
                end
            end
        end

        function p = initMatFileObject(sys,p)
            % Initialize matfile object
            
            matObj = matfile(p.pFilename,'Writable',true);
            matObj.fs = sys.SampleRate;
            matObj.fc = sys.CenterFrequency;
            matObj.userdata = sys.UserData;
            p.pMatObj = matObj;
        end
        
        function p = checkInputData(sys,X,p) %#ok<INUSL>
            % Check input data, cache attributes:
            
            validateattributes(X,{'numeric'}, ...
                {'nonsparse','2d'},mfilename('class'),'X');
            
            sz = size(X);
            spf = sz(1);
            p.pNumChannels = sz(2);
            p.pSamplesPerFrame = spf;
            
            % Very important! Assume all data will be complex in general.
            % The first sample may be real (say, pure zero) for accidental
            % reasons.  We cannot change to complex on-the-fly, so force it
            % now.
            p.pLikeValue = zeros([1 1],'like',complex(X(1,1)));
        end

        function p = checkFileNameAndDeleteFile(sys,p)
            % Check filename specified in public property.
            % Delete file if it exists.
            % Store name with extension.
            %
            % Return updated parameter struct.
            
            % If no extension, append .bb to filename.
            % This means we never create an extension-less file name.
            fname = sys.Filename;
            [~,~,e] = fileparts(fname);
            if isempty(e)
                fname = [fname '.bb'];
            end
            % Remove existing baseband file:
            if exist(fname,'file')
                delete(fname);
            end
            
            p.pFilename = fname;
            p.pFullName = fname; % which(fname);
        end
        
        function resetCache(sys)
            
            % Reset cache data and indices, so the next step-call will
            % write the first sample to the cache and to the file.
            p = sys.pInternal;
            
            % Do NOT reset pNeedFlush to false.
            %
            % If a flush was needed, then it's still needed, and the
            % newly-reset dataset needs to be written to the file (possibly
            % overwriting the file as well).
            
            % Reset cache:
            p.pLastSampleWritten = 0; % could decrease to earlier times
            p.pMaxSampleWritten = 0;  % highest time index ever recorded
            p.pTotalSampleCount = 0;
            
            % Set first and last time sample represented by the cache.
            %
            % This sample range doesn't indicate the state of writing data
            % into the cache!  It only defines the range of time samples
            % represented in the cache.  The range of samples changes as
            % the cache is used for later time steps.
            %
            Ns = p.pSamplesInCache;
            p.pCacheSampleRange = [1 Ns]; % [firstCacheSample, lastCacheSample]
            
            % NOTE: Complexity gets set via the 'like' option.
            p.pCacheSamples = zeros(Ns,p.pNumChannels,'like',p.pLikeValue);
            
            sys.pInternal = p;
        end

        function resetCacheAndFile(sys)
            
            resetCache(sys);
            
            % Remove data previously written to the MAT file.
            %
            % Use correct data type.
            %
            % Write a scalar, not an empty - the reason is that empty does
            % not retain the complex attribute, and when written by the
            % MatFile object, it writes as a real.  Subsequent
            % complex-valued writes will throw an error due to attribute
            % mismatch.  A scalar retains the attribute.
            %
            % This is only an inefficiency and NOT a corruption in general,
            % as this scalar value will be OVERWRITTEN on the next write,
            % since an absolute index is used and the first write will
            % cause that to be 1 (i.e., pLastSampleWritten is now 0).
            %
            % HOWEVER!  Closing the file after reset would show a scalar
            % zero, hence the close(sys) code checks pLastSampleWritten and
            % cleans up this issue.
            %
            % NOTE: Complexity gets set via the 'like' option.
            p = sys.pInternal;
            p.pMatObj.data = zeros(1,p.pNumChannels,'like',p.pLikeValue);
            sys.pInternal = p;
        end
        
        function stepMaxShortData(sys,X)
            % Finite file size.
            % Write only sys.MaxSamples number of samples to file.
            % Frame size <= MaxSamples.
            % This is the typical finite-size situation.
            %
            % Remap indices:
            % Nmax = 10, F = 5
            % [1 2 3 4 5 6 7 8 9 10] <- MaxSize = 10 samples
            %        L               <- Last write position
            %          1 2 3 4 5     <- New samples written
            % [1 2 3 4 5 6 7 8 9 10]
            %            L           <- Last write position
            %  5           1 2 3 4   <- New samples written
            
            % All samples get written.
            p = sys.pInternal;
            F = size(X,1);
            Nmax = sys.MaxSamples;
            i1 = 1+rem(p.pLastSampleWritten,Nmax); % add 1 to last write
            i2 = 1+rem(i1+F-2,Nmax);
            if i1 < i2  % one write will do
                writeData(sys,X,i1);
            else        % two writes, data wraps around finite buffer
                P = Nmax-i1+1;
                writeData(sys,X(1:P,:),i1);
                writeData(sys,X(P+1:end,:),1);
            end
        end
        
        function stepMaxLongData(sys,X)
            % Finite file size.
            % Write only sys.MaxSamples number of samples to file.
            % Frame size > MaxSamples.
            % This is a very unusual case.
            %
            % Remap indices:
            % Only Nmax samples get written.
            p = sys.pInternal;
            F = size(X,1);
            Nmax = sys.MaxSamples;
            R = rem(F,Nmax);
            % o1: sample index into data (output) for writing 2nd segment
            %     of input frame
            o1 = 1+rem(p.pLastSampleWritten+R,Nmax);
            if o1==1
                % All samples, in order
                writeData(sys,X(F-Nmax+1:F,:),1);
            else
                % Two segments to write into buffer:
                d = Nmax-o1+1;
                i2 = F-(Nmax-d)+1;
                i1 = i2 - d;
                writeData(sys,X(i2:F,:),1);
                writeData(sys,X(i1:i2-1,:),o1);
            end
        end
        
        function writeData(sys,X,i1)
            % Write data samples.
            %  X: frame data, rows=time, cols=chans.
            %  i1: absolute time index (optional input, 1-based integer.
            %      Assumes next sample to write immediately follows last if
            %      omitted.
            %
            % Assumes forward traversal of data index requests (ie,
            % increaing i1 on each successive call)
            %
            % Writes data to cache with priority, then to file as needed.
            
            p = sys.pInternal;
            if nargin < 3
                % No input i1; assume the next sample to write is
                % immediately adjacent to the last data sample:
                i1 = p.pLastSampleWritten + 1;
            end
            
            % Write cache samples and clear out cache.
            %
            % writeData(sys,X,i1) writes X(i1:i2) to the
            % cache, if it overlaps - and if it does overlap, it's only
            % partial.  That is written to file, the cache is reset, and
            % the remaining data in X(i1:i2) that was not written to the
            % last cache page is written to the newly reset cache page.
            %
            % Index i1 to i2 are being requested to write
            % Index i2 could not be written to cache (insufficient room)
            % Write cache content to file since it is now full
            % Cache iN+1:i2 in first part of cache,
            %     plus empty cache space for future writes
            %
            % Note that i2 is usually i1+FrameSize-1, but may be less than
            % that on a final partial-frame write-request.
            
            % Determine input data time-stamp range, i1:i2
            F = size(X,1);
            i2 = i1 + F-1;
            p.pTotalSampleCount = p.pTotalSampleCount + F;
            p.pNeedFlush = p.pNeedFlush || F>0;

            % CSR: time range represented by cache
            %   [firstCacheSample lastCacheSample]
            csr = p.pCacheSampleRange;
            fcs = csr(1); % firstCacheSample
            lcs = csr(2); % lastCacheSample
            
            % Determine if any samples can fit in cache
            %
            % Allow time-gap between data sequences
            %   - all samples in X are consecutive, but
            %   - i1 might not be pLastSampleWritten+1
            %
            if i1 >= fcs && i2 < lcs
                % Cache hit: ALL new data samples fit in cache page.
                %
                % If it's an "exact fill" (i2==lcs), we don't enter here
                % and instead do the longer job of cache write below.
                
                p.pCacheSamples(i1+1-fcs:i2+1-fcs,:) = X;
                p.pLastSampleWritten = i2;
                p.pMaxSampleWritten = max(p.pMaxSampleWritten,i2);
                
                % For speed, don't do anything more to cache.
                sys.pInternal = p;
                return % EARLY EXIT
            end
            
            if i1 >= fcs && i1 <= lcs
                % Write first few samples to last few remaining spots in
                % cache.
                %
                % Bump write-request to one sample past end of cache.
                p.pCacheSamples(i1+1-fcs : end,:) = X(1:1+lcs-i1,:);
                iX = 1+lcs-i1+1; % next relative sample in X that needs to be copied
                i1 = 1+lcs;      % next absolute time index corresponding to iX
            else
                % No new-data timestamps are located in cache.
                % Cache needs to be written so a new cache page can begin.
                %
                iX = 1; % next relative sample in X that needs to be written
                % i1: Leave as-is, corresponds to first sample in X (iX=1)
            end
            
            % Write full cache to matfile
            Nc = p.pSamplesInCache;  % # samples per frame
            if p.pFastWrite
                % For more efficient writes, use:
                %    matlab.internal.language.partialLoad(src,data,var)
                % where
                %   src: fully-qualified MAT file name
                %   var: vsubset with the fields:
                %      .Name = 'data'
                %      .Subsets.Type = '()'
                %              .Index = {[sample1 sampleStep sample2], ...
                %                        [chan1 chanStep chan2]}
                vsubset = p.pPartialSaveStruct;
                vsubset.Subsets.Index{1} = [fcs lcs];
                matlab.internal.language.partialSave( ...
                    p.pFullName, p.pCacheSamples,vsubset);
            else
                % Always specify 2nd dim explicitly on LHS - don't just use
                % a colon. The reason is subtle: if the MAT file gets
                % deleted while the object remains in the workspace, then a
                % subsequent step() is performed, pMatObj will try to do
                % the right thing.  But a "colon" for the 2nd dim is
                % underspecified for a new file.  Supplying explicit
                % indices will allow pMatObj to create a new MAT file
                % mid-stream.
                %
                p.pMatObj.data(fcs:lcs, 1:p.pNumChannels) = p.pCacheSamples;
            end
            
            % Move cache page to another sample region
            %
            % - we need to determine if we're revisiting a sample range
            %   that we've already written to
            % - this happens for finite MaxSamples, where a circular buffer
            %   forces us to revisit early samples periodically
            % - if we're revisiting an old range, we need to RELOAD those
            %   samples (in case we're skipping time and some old samples
            %   will still be retained)
            %     xxx this seems unlikely, possibly undesired
            %         if we revisit old, it's only due to time-wrap
            %
            % - if we're visiting a new range, we simply load ZEROS
            %
            csr = [i1 i1+Nc-1];
            p.pCacheSampleRange = csr; % [i1 i1+Ns-1];
            
            % NOTE: Complexity gets set via the 'like' option.
            p.pCacheSamples = zeros(Nc,p.pNumChannels,'like',p.pLikeValue);
            
            % Overwrite start of cache with remaining input data
            Nx = F-iX+1; % # samples in frame of X remaining uncopied
            if Nx > 0 % otherwise can allow out-of-range iX
                p.pCacheSamples(1:Nx,:) = X(iX:end,:);
            end
            
            % Update indices for next write
            lsw = i1+Nx-1;
            p.pLastSampleWritten = lsw;
            p.pMaxSampleWritten = max(p.pMaxSampleWritten,lsw);
            
            sys.pInternal = p;
        end
        
        function flushCache(sys)
            % Flush cache to MAT file and reset cache.
            %
            % None, or only a portion, of the cache page may be filled.
            % Only write to file data actually written to the cache.
            
            p = sys.pInternal;
            if isempty(p) || ~p.pNeedFlush
                % setupImpl never ran prior to flush
                return
            end
            %idx = p.pLastSampleWritten;
            msw = p.pMaxSampleWritten;
            if msw==0
                % This was a "reset" condition without any prior step calls
                % (writes).
                %
                % On reset, we write a zero to the MAT file, due to
                % complexity-attribute issues, and that value gets
                % overwritten in subsequent write operations.  If no writes
                % occurred since reset, we need to clean up that zero now.
                %
                % NOTE: A MAT file is created in this case.  Since some
                % data was written, and then a reset occurred, the user
                % expects a file with EMPTY data in it.
                %
                % This is UNLIKE the situation of having instantiated this
                % object without any step operations, then closing.  In
                % that case, a file is NOT generated - for that case, we
                % return early for an empty pInternal data struct (as
                % detected above).
                %
                % NOTE: Complexity gets set via the 'like' option.
                p.pMatObj.data = zeros(0,p.pNumChannels,'like',p.pLikeValue);
                
            else
                % We've used the cache, so we now check the current cache
                % page -- but just the current cache page.
                %
                % If pLastSampleWritten is anywhere in this page, we flush
                % the page.  That suffices for ALL samples written with
                % this page --- if any sample was written that was NOT in
                % this page, we'd have a different page!  Bottom line: we
                % can just check pLastSampleWritten.
                %
                lsw = p.pLastSampleWritten;
                csr = p.pCacheSampleRange;
                fcs = csr(1); % firstCacheSample
                lcs = csr(2); % lastCacheSample
                
                if lsw>=fcs && lsw<=lcs 
                    % Don't write more of the cache than was written into
                    % it -- otherwise we will copy bad values.  Especially
                    % important when finite size wrapping is used.
                    %
                    i1 = fcs; % first sample written in cache, absolute time index
                    i2 = lsw; % last sample written, absolute time index
                    %i2 = min(lcs,msw);
                    N = i2-i1+1;
                    
                    % Write cache samples to matfile
                    if p.pFastWrite
                        % For more efficient writes, use:
                        %    matlab.internal.language.partialLoad(src,data,var)
                        % where
                        %   src: fully-qualified MAT file name
                        %   var: vsubset with the fields:
                        %      .Name = 'data'
                        %      .Subsets.Type = '()'
                        %              .Index = {[sample1 sampleStep sample2], ...
                        %                        [chan1 chanStep chan2]}
                        vsubset = p.pPartialSaveStruct;
                        vsubset.Subsets.Index{1} = [i1 i2];
                        matlab.internal.language.partialSave(p.pFullName, ...
                            p.pCacheSamples(1:N,:),vsubset);
                    else
                        p.pMatObj.data(i1:i2,:) = p.pCacheSamples(1:N,:);
                    end
                end
            end
        end
        
        function linearizeCircularData(sys)
            % File may contain "circular buffer" data when a long data
            % stream is written to a file where a finite MaxSample value is
            % specified.  In that case, data may be in two contiguous
            % chunks and is reordered to make time continuous from lowest
            % to highest timestamp.
            %
            % There is typically a lot of data to move, and files are used
            % instead of memory to carry out the reordering.
            
            % If maximum timestamp of data written is < MaxSamples, the
            % buffer content never had a chance to wrap, and it's all in
            % order - there's nothing to do.
            %
            % Also, if step never ran, pInternal never gets set, and it
            % remains empty.  For this case, we can also return early.
            Nmax = sys.MaxSamples;
            p = sys.pInternal;
            if isinf(Nmax) || isempty(p) || ~p.pNeedFlush
                return
            end
            Wmax = p.pMaxSampleWritten;
            if Wmax < Nmax
                return
            end
            
            % We copy the two out-of-order segments of the current file to
            % be in-order in a new MAT file, then replace the original file
            % with this new MAT file.
            %
            % Create a new MAT file:
            fname_tmp = [tempname '.mat'];
            mNew = matfile(fname_tmp,'writable',true);
            
            % Get current MAT file (with circular buffer segments):
            p = sys.pInternal;
            mCurr = p.pMatObj;
            
            %  Stream the second contiguous segment (samples N+1:end) to
            %  the start of the new MAT file.  This is "early timestamp"
            %  segment.
            Nt = min(sys.MaxSamples,p.pTotalSampleCount);
            brkpt = p.pLastSampleWritten; % last sample of segment1
            frameSize = p.pSamplesInCache; % frame size used for streaming
            Nc = p.pNumChannels;
            src1 = brkpt+1; % oldest sample in circular buffer
            src2 = src1;
            outOffset = 0;
            while src2 < Nt
                src2 = min(src1+frameSize-1,Nt);
                N = src2-src1+1;
                dst1 = outOffset+1;
                dst2 = outOffset+N;
                mNew.data(dst1:dst2,1:Nc) = mCurr.data(src1:src2,:);
                outOffset = outOffset + N;
                src1 = src2+1; % next starting sample to copy
            end
            
            %  Stream the first contiguous segment (samples 1:N) to
            %  the end of the new MAT file.  This is "late timestamp"
            %  segment.
            %Nt = p.pTotalSampleCount; % total
            brkpt = p.pLastSampleWritten; % last sample of segment1
            frameSize = p.pSamplesInCache;
            Nc = p.pNumChannels;
            src1 = 1;
            src2 = 1;
            while src2 < brkpt
                src2 = min(src1+frameSize-1,brkpt);
                dst1 = outOffset+src1;
                dst2 = outOffset+src2;
                mNew.data(dst1:dst2,1:Nc) = mCurr.data(src1:src2,:);
                src1 = src2+1; % next starting sample to copy
            end
            
            % Move new MAT file to replace old MAT file:
            src = mNew.Properties.Source;
            dst = sys.pInternal.pFullName;
            [success,msg,msgid] = movefile(src,dst,'f');
            if ~success
                error(msgid,msg);
            end
        end
    end
end

function bps = bytesPerSample(val) %#ok<INUSD>
% Return # bytes per element.
% Accounts for complex values.
s = whos('val');
bps = s.bytes;
end
