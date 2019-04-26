classdef BasebandReader < matlab.System & matlab.system.mixin.FiniteSource
    %BasebandReader Return baseband data from a compressed file.
    %   S = BasebandReader returns a System object S that reads complex-
    %   valued data from a compressed data file.  One frame of data is
    %   returned at a time.
    %
    %   S = BasebandReader('PropertyName',PropertyValue, ...) returns a
    %   System object with named properties set to specified values.
    %
    %   S = BasebandReader(FNAME,SPF,'PropertyName',PropertyValue,...)
    %   returns a BasebandReader object where 'Filename' is set to FNAME
    %   and 'SamplesPerFrame' is set to SPF. If a file extension is omitted
    %   from FNAME, '.bb' is assumed.  These values can also be set using
    %   property name and value pairs.
    %
    %   A baseband data file is a specific type of MAT file written by
    %   BasebandWriter. Baseband data is a matrix of complex-valued (I/Q)
    %   samples of a signal that have been down-converted from a center
    %   frequency of FC Hz to 0 Hz ("baseband"), and sampled at a rate of
    %   FS Hz. The matrix contains consecutive time samples down each
    %   column, and one column per data channel. FC and FS are recorded
    %   when a baseband data file is created, and may be accessed while the
    %   file is being read from the read-only properties CenterFrequency
    %   and SampleRate. Additional meta-data stored when the baseband file
    %   was written by BasebandWriter may be obtained from the UserData
    %   property.
    %
    %   Samples are returned from the BasebandReader object when the step
    %   method is called:
    %      BB = step(S);
    %   BB is a SamplesPerFrame-by-NumChannels matrix of complex values. To
    %   identify if the last sample of the file has been returned in BB, an
    %   additional return argument EOD may be requested:
    %     [BB,EOD] = step(S);
    %
    %   After the last data value is read from the file, zeros are
    %   returned. Call reset(BB) to rewind to the start of the file at any
    %   time.
    %
    %   isDone(S) returns true once the last data sample has been
    %   returned.  In addition, SamplesRead is a property that returns
    %   the number of samples read from the baseband file.
    %
    %   See also BasebandWriter.

    % Copyright 2015 The MathWorks, Inc.
    
    properties (Nontunable)
        % Name of baseband data file.  If a file extension is omitted,
        % '.mat' is assumed.
        Filename = 'data'
    end
    
    properties (SetAccess=private)
        % Center frequency of RF signal prior to downconversion to
        % baseband, in Hz.
        CenterFrequency = 0
        
        % Sample rate of complex baseband samples, in Hz.
        SampleRate = 0
        
        % Number of columns in data matrix
        NumChannels = 0
    end

    properties (Nontunable)
        % Number of samples of baseband data to return, per channel, from
        % each call to the step function.
        SamplesPerFrame = 64
    end
    
    properties (SetAccess=private,Dependent)
        % Index of last sample read from file.
        SamplesRead
    end
    
    properties (Nontunable,Dependent,Hidden)
        % Number of frames to read from data file at a time. Larger values
        % increase performance at the expense of increased memory usage.
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
    
    properties (SetAccess=private)
        % User-specified data stored with the object.
        UserData
    end
    
    properties (Access=private)
        % Struct of internal properties, used for efficiency.  It is faster
        % to access struct fields versus individual object properties.
        %   .pFastRead
        %   .pMatObj
        %   .pLikeValue
        %   .pSize
        %   .pFilename
        %   .pFullName          (empty unless pFastRead=true)
        %   .pPartialLoadStruct (empty unless pFastRead=true)
        %   .pSamplesPerFrame
        %   .pCacheSampleRange
        %   .pCacheSamples
        %   .pFramesInCache
        %   .pWhosInfo
        %   .pLastSampleRead
        pInternal
        
        % Used to hold manual value of CacheSize.  Needed so automatic
        % value computation doesn't trigger 'manual' mode.
        pCacheSize = 1
    end
    
    methods
        function sys = BasebandReader(varargin)
            setProperties(sys,nargin,varargin{:}, ...
                'Filename','SamplesPerFrame');
            
            % In case default file is going to be used
            init(sys); % side-effect
        end
        
        function set.Filename(sys,val)
            if ~ischar(val)
                error('Filename must be a string.');
            end
            sys.Filename = val;
            init(sys); % side-effect
        end
        
        function set.SamplesPerFrame(sys,val)
            validateattributes(val,{'numeric'}, ...
                {'scalar','real','positive','integer','finite'}, ...
                mfilename('class'),'SamplesPerFrame');
            sys.SamplesPerFrame = val;
            init(sys); % side-effect
        end
        
        function y = get.CacheSize(sys)
            y = sys.pCacheSize;
        end
        
        function set.CacheSize(sys,val)
            validateattributes(val,{'numeric'}, ...
                {'scalar','real','positive','integer','finite'}, ...
                mfilename('class'),'CacheSize');
            
            sys.pCacheSize = val;
            
            % This side-effect triggers a call to init:
            sys.CacheSizeMode = 'manual';
        end
        
        function set.CacheSizeMode(sys,val)
            sys.CacheSizeMode = validatestring(val, ...
                {'auto','manual'}, mfilename('class'),'CacheSizeMode');
            
            init(sys); % side-effect
        end
        
        function y = get.SamplesRead(sys)
            p = sys.pInternal;
            if isempty(p)
                y = 0;
            else
                y = sys.pInternal.pLastSampleRead;
            end
        end
    end
    
    methods (Access=protected)
        function setupImpl(sys)
            init(sys,true);
        end
        
        function y = isDoneImpl(sys)
            p = sys.pInternal;
            y = p.pLastSampleRead >= p.pSize(1);
        end
        
        function resetImpl(sys)
            % Reset so the next step will return the first sample from the
            % file.
            p = sys.pInternal;
            p.pLastSampleRead   = 0;
            p.pCacheSampleRange = [0 0];
            p.pCacheSamples     = []; % data samples in cache
            sys.pInternal = p;
        end
        
        function ret = infoImpl(sys)
            % Format data struct for "info" display
            
            p = sys.pInternal;
            if isempty(p)
                s.ComplexData = [];
                s.DataClass = '';
                s.NumChannels = [];
                s.NumFrames = [];
                s.NumSamples = [];
                s.BytesInCache = [];
                s.FramesInCache = [];
            else
                w = p.pWhosInfo;
                s.ComplexData = w.complex;
                s.DataClass = w.class;
                sz = w.size;  % = size(matObj,'data');
                s.NumChannels = sz(2);
                s.NumFrames = ceil(sz(1)/sys.SamplesPerFrame);
                s.NumSamples = sz(1);
                s.BytesInCache = p.pFramesInCache * p.pSamplesPerFrame ...
                    * bytesPerSample(p.pLikeValue);
                s.FramesInCache = p.pFramesInCache;
            end
            
            if nargout > 0
                ret = s;
            else
                disp(s);
            end
        end
        
        function [y,isDone] = stepImpl(sys)
            % Read next frame of samples from baseband data file.
            %
            % isDone is true if samples past end of file are returned.
            
            p = sys.pInternal;
            i1 = 1 + p.pLastSampleRead;
            sz = p.pSize;
            if i1 > sz(1)
                % Past end of data.  All samples are zero padding.
                %
                % No need to update pLastSampleRead
                y = zeros([p.pSamplesPerFrame p.pSize(2)],'like',p.pLikeValue);
                isDone = true;
            else
                i2 = i1-1 + p.pSamplesPerFrame;
                isDone = i2 > sz(1);
                if isDone
                    % Past end of data. Some samples are zero padding.
                    % - use getData for samples remaining in file/cache
                    % - append zeros to fill a full frame
                    % - update pLastSampleRead
                    y = [getData(sys,i1,sz(1)); ...
                        zeros([i2-sz(1) sz(2)],'like',p.pLikeValue)];
                else
                    % Full data frame
                    % - no appending of zeros needed
                    % - update pLastSampleRead
                    y = getData(sys,i1,i2);
                end
            end
        end
    end
    
    methods (Access=private)
        function p = clearInternals(sys,p)
            
            sys.CenterFrequency = 0;
            sys.SampleRate = 0;
            sys.NumChannels = 0;
            sys.UserData = [];
            
            % sys.SamplesRead is Dependent, reset this instead:
            p.pLastSampleRead = 0;
        end
        
        function init(sys,throwError)
            % Initialization performed whenever a non-tunable parameter
            % changes value.
            
            if nargin < 2
                throwError = false;
            end
            
            % Choose MAT-file API:
            %   false: public API, slower file-reads
            %    true: internal-only API, faster file-reads
            p.pFastRead = false;
            
            [p,isValid] = checkFileName(sys,p,throwError);
            if isValid
                p = initMatFileObject(sys,p);
                p = initCacheSize(sys,p);
                p = initFastRead(sys,p);
            else
                p = clearInternals(sys,p);
            end
            % Use struct of private properties for access speed
            sys.pInternal = p;
        end
        
        function y = getData(sys,i1,i2)
            % getData(sys,i1,i2) retrieves data samples from index i1
            % through i2 inclusive.
            %
            % Assumes forward-traversal of data index requests, i.e., that
            % i1 and i2 are increasing since the last call to getData.
            %
            % Pulls from cache with priority, then from file as needed.
            % Updates pLastSampleRead.
            
            % Lookup samples in cache
            p = sys.pInternal;
            csr = p.pCacheSampleRange; % [firstCacheSample lastCacheSample]
            fcs = csr(1);              % firstCacheSample
            if i1 >= fcs && i2 <= csr(2)
                % Cache hit
                y = p.pCacheSamples(i1+1-fcs:i2+1-fcs,:);
            else
                % Cache miss (total or partial)
                y = fetchAndFillCache(sys,p,csr,fcs,i1,i2);
            end
            sys.pInternal.pLastSampleRead = i2;
        end
        
        function y = fetchAndFillCache(sys,p,csr,fcs,i1,i2)
            % Fetch samples and refill cache.
            %
            % Index i1 to i2 are being requested.
            % Index i1 might be in the cache.
            % Index i2 is not in the cache.
            %
            % Fetch i1:i2 for immediate return.
            % Cache i2+1 : i2+1+lookAhead
            %
            % Note that i2 is usually i1+FrameSize-1, but may be less than
            % that on a final partial-frame request.
            
            % Determine if any samples are already in cache
            lcs = csr(2); % lastCacheSample
            if i1 >= fcs && i1 <= lcs
                % Use remaining samples from cache.
                % Bump read-request to one sample past end of cache.
                priorCacheSamples = p.pCacheSamples(i1+1-fcs : end,:);
                i1 = lcs+1;
            else
                % Complete miss:
                %   i1 is before start of cache
                % or
                %   i1 is past end of cache
                priorCacheSamples = [];
            end
            
            % Determine read-request range.
            % Use "absolute" index coordinates.
            % Do not read past end of data!
            r1 = i1;
            r2 = min(i2 + p.pFramesInCache * p.pSamplesPerFrame, ...
                p.pSize(1));
            
            % Read data samples from matfile object:
            if p.pFastRead
                % For more efficient reads, use:
                %    y = matlab.internal.language.partialLoad(src,var,'-mat')
                % where
                %   src: fully-qualified MAT file name
                %   var: struct with the fields:
                %      .Name = 'data'
                %      .Subsets.Type = '()'
                %              .Index = {[sample1 sampleStep sample2], ...
                %                        [chan1 chanStep chan2]}
                %
                vsubset = p.pPartialLoadStruct;
                % Option 1: This requires removal of the private-property
                %           attribute from the following classes:
                %      matlab.internal.language.Subset: .Index
                %      matlab.internal.language.VariableSubset: .Subsets
                %
                vsubset.Subsets.Index{1} = [r1 r2];
                %
                % vsubset.TunableSubsets.TunableIndex{1} = [r1 r2];
                %
                % Option 2: Alternatively, this requires adding a hidden
                %           method to each class:
                %      matlab.internal.language.Subset: updateIndex()
                %      matlab.internal.language.VariableSubset:
                %                                     updateSubsets()
                %
                %index = vsubset.Subsets.Index;
                %index{1} = [r1 r2];
                %vsubset = updateSubsets(vsubset,updateIndex(vsubset.Subsets,index));
                
                % Return the data:
                data = matlab.internal.language.partialLoad(p.pFullName, ...
                    vsubset,'-mat');
            else
                % Slower reads, but supported without modification to the
                % existing Subset and VariableSubset classes:
                data = p.pMatObj.data(r1:r2,:);
            end
            
            % Get subset of data to fulfill current request - i1:i2
            Ncurr = i2-i1+1;
            y = [priorCacheSamples; data(1:Ncurr,:)];
            
            % Retain remainder of data in cache after the return of y.
            %
            % Only retain i2+1:r2, since we're returning i1:i2 now.
            %
            % That will be a full pFramesInCache number of cached frames.
            %
            p.pCacheSampleRange = [i2+1 r2]; % absolute indices
            p.pCacheSamples = data(Ncurr+1:end,:); % a full cache
            sys.pInternal = p;
        end
        
        function [p,isValid] = checkFileName(sys,p,throwError)
            % Check baseband file name
            % Store with extension in pInternal
            
            fname = sys.Filename;
            isValid = exist(fname,'file');
            if ~isValid
                % If no extension, append .bb to file and try again
                [~,~,e] = fileparts(fname);
                if isempty(e)
                    fname = [fname '.bb'];
                    isValid = exist(fname,'file');
                end
            end
            if ~isValid && throwError
                error('comm:baseband:FileNotFound', ...
                    'Baseband file "%s" not found.',sys.Filename);
            end
            p.pFilename = fname;
        end

        function throwIncompatibleFileError(sys,optStr,fname)
            % Throw error indicating how to convert legacy MAT file to a
            % new baseband MAT file.
            %
            % Uses: sys.pInternal.pFilename (if fname not passed)
            
            if nargin < 3
                fname = sys.pInternal.pFilename;
            end
            if nargin < 2
                optStr = '';
                optFmt = '%s';
            end
            if ~isempty(optStr)
                optFmt = '(%s)\n';
            end
            error('comm:baseband:InvalidFile', ...
                ['File ''%s'' was not created by BasebandWriter.\n' ...
                optFmt ...
                '\nConvert it to a compatible file using:\n' ...
                'BasebandWriter.convertFile(''%s'',''VarName'',''NewFileName'',Fs,Fc);'], ...
                fname,optStr,fname);
        end
        
        function p = initMatFileObject(sys,p)
            % Create MATFILE object
            % - test variables contained within it
            % - cache variable attributes
            %
            % Uses:
            %       p.pFilename
            % Sets:
            %       p.pMatObj
            %       p.pWhosInfo
            %       p.pSize
            %       p.pLikeValue
            %       sys.SampleRate
            %       sys.CenterFrequency
            
            fname = p.pFilename;
            matObj = matfile(fname,'Writable',false);
            
            if ~matObj.Properties.SupportsPartialAccess
                % Specified file not saved using the -v7.3 option, so it
                % can't be a baseband file:
                optStr = 'File was not saved using the -v7.3 option.';
                throwIncompatibleFileError(sys,optStr,fname);
            end
            p.pMatObj = matObj;
            
            % Remove automatic 'Properties' property - it's not data from
            % the file, it's a matfile service:
            %{
            pnames = properties(matObj);
            sel = strcmpi(pnames,'Properties');
            pnames = pnames(~sel);
            %}
            %
            % Alternative:
            pnames = who(matObj);
            
            % Check for required variable names
            reqd_vars = {'data','fc','fs'};
            vars_present = ismember(reqd_vars,pnames);
            if ~all(vars_present)
                % Build string of missing variables: "'data','fs'"
                missing_cstr = reqd_vars(~vars_present);
                missing_str = sprintf('''%s'',',missing_cstr{:});
                missing_str = missing_str(1:end-1); % remove last comma
                optStr = ['File is missing required variables: ' ...
                    missing_str];
                throwIncompatibleFileError(sys,optStr,fname)
            end
            
            % Get data-variable info
            w = whos(matObj,'data');
            p.pWhosInfo = w;
            
            if w.sparse
                optStr = 'File contains sparse data which is unsupported.';
                throwIncompatibleFileError(sys,optStr,fname)
            end
            
            % Cache required info:
            sys.SampleRate = matObj.fs;
            sys.CenterFrequency = matObj.fc;
            
            % Size of file contents, [total # samples, # channels]
            sz = w.size;  % = size(matObj,'data');
            p.pSize = sz;
            sys.NumChannels = sz(2);
            
            % Get optional userdata info from MAT file:
            if any(strcmpi('userdata',pnames))
                sys.UserData = matObj.userdata;
            else
                sys.UserData = [];
            end
            
            % Copy a scalar value from data, for attributes.
            %
            % NOTE: Partial access will trigger checking of the MAT
            % file format, throwing an error if it is not -v7.3:
            scalarVal = matObj.data(1,1);
            p.pLikeValue = zeros([1 1],'like',scalarVal);
        end
        
        function p = initCacheSize(sys,p)
            
            % Properties like SamplesRead are visible after construction,
            % and we want it to show 0 and not empty:
            p.pLastSampleRead = 0;
            
            % Copy some read-only property values for faster access
            spf = sys.SamplesPerFrame;
            p.pSamplesPerFrame = spf;
            
            % Calculate # frames in file, rounded up so last sample in file
            % will be in cache.  Therefore last frame may be only partially
            % filled.
            sz = p.pSize;
            Nf_file = ceil(sz(1)/spf);
            
            if strcmpi(sys.CacheSizeMode,'auto')
                % Automatic cache size (all in units of a frame!)
                %    = max(1,min(256kB,FileSize))
                bpf = bytesPerSample(p.pLikeValue) * spf;
                Nf_256k = floor(256*1024 / bpf);
                sys.pCacheSize = max(1,min(Nf_256k,Nf_file));
            end
            
            % Limit cache size to max # frames in file
            p.pFramesInCache = min(Nf_file,sys.pCacheSize);
            
            % Setup data read cache:
            p.pCacheSampleRange = [0 0];
            p.pCacheSamples = [];
        end
        
        function p = initFastRead(sys,p) %#ok<INUSL>
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
            if p.pFastRead
                sz = p.pSize;
                subset = matlab.internal.language.Subset('()',{[1 1],[1 sz(2)]});
                p.pFullName = which(p.pFilename);
                p.pPartialLoadStruct = ...
                    matlab.internal.language.VariableSubset('data',subset);
            else
                p.pFullName = '';
                p.pPartialLoadStruct = [];
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
