% SHOW Prints the full contents of a multi-element or array structure
%
% Description: Prints the full contents of a structure variable to the Command
%              Window or to a file. Ideal for showing contents of a nested
%              structure or a structure array.
%
% Motivation:
%     Often it takes multiple commands to explore the contents of a structure,
%     especially if it is a structure array or has nested structure elements.
%     For example, take this relatively simple case:
%       >> s = struct('data',{struct('value',{0,'a',pi}),'str',magic(3),1-2i})
%     This structure is displayed by MATLAB as:
%       s =
%       1x4 struct array with fields:
%           data
%     Not much information there, other than the top level size and field(s).
%     To explore the actual contents, you would need to type:
%       >> s(1)
%       >> s(1).data(1)
%       >> s(1).data(2)
%       >> s(1).data(3)
%       >> s(2)
%       >> s(3)
%       >> s(4)
%     With SHOW, all that effort is done for you in one simple command:
%       >> show(s)
%     Also, the printed output is more readable and concise, with field names
%     tab-indented to represent nested structure elements.
%
% Inputs:
%     VAR  - (MxN struct) variable to be shown
%     FID  - (optional) file identifier to print contents to a file
%
% Outputs: None.
%
% Usage:
%     show(var)
%       -or-
%     fid = fopen('file.txt','w'); show(var,fid); fclose(fid);
%
% Example:
%     % Show full contents of the image info (needs Image Processing Toolbox)
%     structInfo = imfinfo('cloudCombined.jpg');
%     show(structInfo);
%
% Example:
%     clc
%     myStruct = struct('nestedStruct',struct('version','1.0','size',[800 600], ...
%         'path','/home/user/matlab','file','data.mat', 'matrixData',magic(3), ...
%         'img',uint8(rand(600,800,3)),'cellField',{{'a','b','c'}}, ...
%         'shortStringField','A string with 100 characters or less is printed', ...
%         'longStringField',repmat('Only size/type shown for long strings. ',1,3), ...
%         'rowVectorData',rand(1,3),'columnVectorData',sort(rand(5,1))), ...
%       'structArray',struct('data',{[1+4i,3-2i],randperm(10),struct('val',{'a',2,pi})}))
%     fprintf('\nCompare default MATLAB "disp" above versus "show" below:\n\n');
%     show(myStruct);
%
% Example:
%     % Print contents to a file
%     [uv,sv] = memory;
%     memoryInfo = struct('UserView',uv,'SystemView',sv);
%     fid = fopen('memory.txt','w');
%     show(memoryInfo,fid);
%     fclose(fid);
%
% Example:
%     % Print contents to the screen
%     [uv,sv] = memory;
%     memoryInfo = struct('UserView',uv,'SystemView',sv);
%     show(memoryInfo);
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.0
% Release Date: 10/5/15
%
function show(var,fid,nTabs)
    
    
    % Set default print stream
    if nargin < 2
        fid = 1;
    end
    
    % Set default number of tabs
    if nargin < 3
        nTabs = 1;
    end
    
    % Identify variable information
    nVals    = numel(var);
    tabStr   = repmat('\t',1,nTabs);
    varSize  = size(var);
    varDims  = length(varSize);
    varClass = class(var);
    typeStr  = sprintf(['%d' repmat('x%d',1,varDims-1) ' %s'],varSize,varClass);
    
    % Print contents based on variable type
    switch varClass
        case 'struct'
            % Show struct contents
            fieldNames = fieldnames(var);
            nFields    = length(fieldNames);
            nameLength = cellfun(@length,fieldNames);
            maxLength  = max(nameLength);
            numSpaces  = maxLength - nameLength;
            % Handle struct arrays
            for i = 1:nVals
                % Show input variable name
                if (nTabs == 1)
                    inputName = inputname(1);
                    if ~isempty(inputName)
                        if (nVals == 1)
                            fprintf(fid,'%s:\n',inputName);
                        else
                            fprintf(fid,'%s(%d):\n',inputName,i);
                        end
                    end
                end
                % Show contents for all struct fields
                for j = 1:nFields
                    thisField = fieldNames{j};
                    thisValue = var(i).(thisField);
                    thisSpace = repmat(' ',1,numSpaces(j));
                    if isstruct(thisValue)
                        nValues = length(thisValue);
                        if (nValues == 1)
                            % Show single-element structure contents
                            fprintf(fid,[tabStr '%s:\n'],thisField);
                            show(thisValue,fid,nTabs+1);
                        else
                            % Handle nested structure arrays
                            for k = 1:nValues
                                fprintf(fid,[tabStr '%s(%d):\n'],thisField,k);
                                show(thisValue(k),fid,nTabs+1);
                            end
                        end
                    else
                        % Show non-structure contents
                        fprintf(fid,[tabStr '%s: ' thisSpace],thisField);
                        show(thisValue,fid,nTabs+1);
                    end
                end
            end
        case 'cell'
            % Show cell contents
            if ~nVals
                fprintf(fid,'{}\n');
            else
                fprintf(fid,'{%s}\n',typeStr);
            end
        case 'char'
            % Show string contents
            if (varSize(1) <= 1) && (nVals <= 100)
                fprintf(fid,'''%s''\n',var);
            else
                fprintf(fid,'[%s]\n',typeStr);
            end
        case {'single','double','logical','int8','uint8', ...
                'int16','uint16','int32','uint32','int64','uint64'}
            % Show numeric contents
            falseTrue = {'false','true'};
            cs = {'+','-'};
            if ~nVals
                fprintf(fid,'[]\n');
            elseif issparse(var)
                typeStr = sprintf(['%d' repmat('x%d',1,varDims-1) ' sparse %s'],varSize,varClass);
                fprintf(fid,'[%s]\n',typeStr);
            elseif (nVals == 1)
                if islogical(var)
                    fprintf(fid,'%s\n',falseTrue{var+1});
                elseif isreal(var)
                    fprintf(fid,'%1.15g\n',var);
                else
                    if imag(var) < 0
                        fprintf(fid,'%1.15g-%1.15gi\n',real(var),abs(imag(var)));
                    else
                        fprintf(fid,'%1.15g+%1.15gi\n',real(var),imag(var));
                    end
                end
            elseif (varSize(1) == 1) && (nVals <= 10)
                if islogical(var)
                    fprintf(fid,['[%s' repmat(' %s',1,nVals-1) ']\n'],falseTrue{var+1});
                elseif isreal(var)
                    fprintf(fid,['[%1.15g' repmat(' %1.15g',1,nVals-1) ']\n'],var);
                else
                    isNeg = imag(var) < 0;
                    fprintf(fid,'[%1.15g%s%1.15gi',real(var(1)),cs{isNeg(1)+1},abs(imag(var(1))));
                    for k = 2:nVals
                        fprintf(fid,', %1.15g%s%1.15gi',real(var(k)),cs{isNeg(k)+1},abs(imag(var(k))));
                    end
                    fprintf(fid,']\n');
                end
            elseif isreal(var)
                fprintf(fid,'[%s]\n',typeStr);
            else
                typeStr = sprintf(['%d' repmat('x%d',1,varDims-1) ' complex %s'],varSize,varClass);
                fprintf(fid,'[%s]\n',typeStr);
            end
        otherwise
            % Show other content types
            fprintf(fid,'[%s]\n',typeStr);
    end
    
end

