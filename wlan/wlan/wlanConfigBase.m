classdef wlanConfigBase
%wlanConfigBase Base object for configuration objects

% Copyright 2015 The MathWorks, Inc.

%#codegen
%#ok<*AGROW>

methods 
  function v = set(obj, prop)
    v = obj.([prop, '_Values']);
  end       
end

methods (Access = protected)    
  function validateEnumProperties(obj, prop, value)
    options = set(obj, prop);             
    isInvalid = ~any(strcmp(value, options));
    
    if isInvalid
        numOptions = length(options);
        optionStr  = cell(1, numOptions);
        for i = 1 : length(options)-2
            optionStr{i} = ['''', options{i}, '''', ', '];
        end
        optionStr{numOptions-1} = ['''', options{numOptions-1}, ''''];
        optionStr{numOptions}   = [' and ', '''', options{numOptions}, ''''];
        
        coder.internal.errorIf(isInvalid, ...
            'wlan:wlanConfigBase:InvalidEnumValue', ['''', prop, ''''], ...
            [optionStr{1:numOptions}]);
    end
  end
  
  function flag = isInactiveProperty(~, ~)
    flag = false;
  end
end

methods
  function obj = wlanConfigBase(varargin)
    coder.internal.errorIf((mod(nargin,2) ~= 0), ...
        'wlan:wlanConfigBase:InvalidPVPairs');

    for i = 1:2:nargin
        obj.(varargin{i}) = varargin{i+1};
    end
  end
  
  function disp(obj)  
    objName = class(obj);
    
    if usejava('swing') && desktop('-inuse') && feature('hotlinks');
        fprintf(['  <a href="matlab:help(', objName, ')">', objName, '</a>', ' with properties:\n\n']);
    else
        fprintf(['  ', class(obj), ' with properties:\n\n']);
    end

    propList = properties(obj);
      
    % Get maximum property name length
    lenBeforeColon = 1;
    activePropList = cell(1, length(propList)); 
    numActiveProp = 0;
    for i = 1:length(propList)
        if ~isInactiveProperty(obj, propList{i})
            numActiveProp = numActiveProp + 1;
            lenBeforeColon = max(lenBeforeColon, length(propList{i}));
            activePropList{numActiveProp} = propList{i};
        end
    end
    % Add an additional tab in front of each line
    lenBeforeColon = lenBeforeColon + 4;
      
    for i = 1:numActiveProp
        prop = activePropList{i};
        thisLineProp = [repelem(' ', lenBeforeColon-length(prop)), prop, ': '];
        if isa(obj.(prop), 'char')                 
            thisLine = [thisLineProp, '''', obj.(prop), '''\n'];
        elseif isa(obj.(prop), 'logical')
            if obj.(prop)
                thisLine = [thisLineProp, 'true\n'];
            else
                thisLine = [thisLineProp, 'false\n'];
            end
        elseif isscalar(obj.(prop))
            thisLine = [thisLineProp, num2str(round(double(obj.(prop))*1e7)/1e7), '\n'];
        elseif isrow(obj.(prop)) && (length(obj.(prop)) <= 4) && all((obj.(prop) == floor(obj.(prop))))
            propVal = obj.(prop);
            rowStr = '[';
            for idx = 1:length(propVal)-1
                rowStr = [rowStr, num2str(floor(propVal(idx))), ' ']; 
            end
            thisLine = [thisLineProp, rowStr, num2str(floor(propVal(end))), ']\n'];
        else % Numeric array 
            propVal = obj.(prop);
            sizeStr = num2str(size(propVal, 1));
            for dimIdx = 2:ndims(propVal)
                sizeStr = [sizeStr, 'x', num2str(size(propVal, dimIdx))]; 
            end
            thisLine = [thisLineProp, '[', sizeStr, ' ', class(propVal), ']\n'];
        end
        fprintf(thisLine);
    end
    fprintf('\n');
  end
end

end

% [EOF]