function [tab,attrout] = xlstblread(filename,sheetname,headerlines,varargin)
%XLSTBLREAD read columwise a mxn table stored in an XLS file (The method uses a robust raw format and applies conversion rules)
%   Syntax: tab = xlstblread(filename,[sheetname],[headerlines],['transpose'],['noheaders'],['emptyisnotnan'])
%           [tab,attributes] = xlstblread(...)
%
%            INPUTS (with the following convention: [] = optional parameter, 'XXX' = keyword)
%          filename: filename (string), list of filenames (cell array of strings), EXPLORE outputs (all EXPLORE formats are accepted)
%                    If missing, EXCEL extensions '.xls' (higher priority) '.xlsx' are automatically added according to available files
%         sheetname: sheet name (string) or list of sheets (cell array of string)
%                    By default, empty cells are replaced by NaN.
%                    Cells, which contain errors, are filled with the error messages and the entire column is assumed to be a cell array.
%                    Non square/rectangle tables are accepted by considering that missing cells are empty.
%                    use 'all' to open all available sheets in the file
%                    use 'first' top open the first sheet only
%                    use 'last' top open the first sheet only
%       headerlines: number of headerlines, where no data are stored
%                    Only the first row is used as variable names (after conversion), the other lines are discarded
%                    When a column title is missing, a title is assigned to this column using 'noheaders' rules
%                    TIP: use [discardedlines headerlines] to discard lines
%       'transpose': keyword (position >3) to transpose the ExcelSheet (row-wise Table instead of column-wise Table)
%       'noheaders': keyword to discard headers as variable names with the following rules
%                    - 'colj' is assigned to the jth column: 'col1',...'coln'
%                    - 'rowj' is applied instead when transpose is used               
%   'emptyisnotnan': prevent empty cells from being replaced by NaNs
%   'notext2values': prevent an attempt the automatic conversion of text into numbers when possible (i.e. when the result of str2double is not NaN)
%'mergeheaderlines': merge header lines
%       'maketable': converts the generated structure into a table and populate attributes as metadata

%
%           OUTPUTS:
%               tab: structure with fields 
%                    tab.(filename).(sheetname).(variablename) = vector
%                    filename, sheetname are assigned only when several files or sheets are used
%                    filename, sheetname, variablename: invalid characters are automatically removed or replaced
%                    variablename is the title of the column (based on the first row)
%                    vector is numerical, when all data in the column are numerical
%                    vector is a cell array otherwise
%                    >For each text cell content, an attempt of conversion in values (double) is performed.
%                    >The conversion is accepted when the result is not a NaN (use 'notext2values' to discard the conversion).).
%                    >Cells, which contains the texts 'NaN' 'Inf','+Inf','-Inf' are considered as numerical (as Matlab does)
%                    >Empty cells are considered as numerical, use 'emptyisnotnan' to discard this feature.
%
%   Examples (from the thesis of Guillaume)
%       tab=xlstblread('c:\data\olivier\e-mail\attach\Chi_database19','complete');  % extract only data from the sheet complete 
%       tab=xlstblread('c:\data\olivier\e-mail\attach\Chi_database19','all');       % extract all sheets
%
%   Dependencies for distribution: EXPLORE, FILEINFO, REPLACEFROMTABLE
%   Character conversion rules can be completed if needed (see variable table in the Definitions section).
%
%   See also: BYKEYWORDS, LOADODS, LOADODSPREFETCH, XLSREAD

% MS 2.0 - 15/01/08 - INRA\Olivier Vitrac - rev. 14/11/15

% Revision History
% 21/01/08 automatic conversion of text into numbers when possible (e.g. after OCR, when mixed types are used)
% 29/01/08 remove '=' from var names and set the max length of var to 16 (VARMAXLEN)
% 29/05/09 fix sheet and file names including figures as first characters
% 18/06/09 fix empty sheets
% 25/08/06 when nargin==1, set sheetname to 'all'
% 04/08/10 add [discardedlines headerlines], remove char(10) and char(13), add 'first' 'last'
% 19/09/15 add mergeheaderlines
% 07/11/15 add attributes, add makedata
% 14/11/15 add '~' to the table of replacement characters ('~' is replaced by '')

% Definitions
kwlist = {'Inf' '+Inf' '-Inf' 'NaN' 'ActiveX VT_ERROR: '}; % values which are replaced 
varprefix = 'col'; % default variable name when no column name is found
varprefixtranspose = 'row'; % variable name when transpose is used
duplicatesuffix = 'dup'; % suffix for duplicated column name
table = {' ' 'é' 'è' 'ê' 'à' 'ù' ':' ',' ';' '.' '-' '~' '+' '*' '\' '/' '°' 'µ' '(' ')' '[' ']' '{' '}' '=' '''' '%' '?' '!' '§' char(13) char(10)
         ''  'e' 'e' 'e' 'a' 'u' ''  ''  ''  ''  ''  '' 'p' 'x' ''  ''  'o' 'u' ''  ''  ''  ''  ''  ''  ''  '_'  'p'  '' ''  ''  ''       ''}'; % add additional character conversion rules if needed
hearderlines_default = 1;
transposeon = false;
allsheetson = false;
firstsheeton= false;
lastsheeton = false;
VARMAXLEN   = 48; % increase this number if required
varlineseparator = '__';
varlineseparatorpat = [regexptranslate('escape',varlineseparator) '$'];

 
%arg check
if nargin<1, error('tab = xlstblread(filename,[sheetname])'), end
if nargin<2, sheetname = ''; end
if nargin<3, headerlines = []; end
if isempty(sheetname), sheetname = 'all'; end
if isstruct(filename), filename = explore(filename,'fullabbreviate'); end
if ~iscell(filename), filename = {filename}; end
if ~iscell(sheetname), sheetname = {sheetname}; end
varargin = lower(varargin);
if strcmpi(sheetname,'all'), allsheetson = true;
elseif strcmpi(sheetname,'first'), firstsheeton = true;
elseif strcmpi(sheetname,'last'), lastsheeton = true;
end
m = length(filename);
n = length(sheetname);
if isempty(headerlines), headerlines = hearderlines_default; end
if ~isnumeric(headerlines), error('headerlines must be an integer>=0 = number of header lines'), end
headerlines = max(0,headerlines);
if ismember('transpose',varargin) || ~any(headerlines), transposeon = true; varprefix = varprefixtranspose; end
nohearderson  = ismember('noheaders',varargin);
emptyisnotnan = ismember('emptyisnotnan',varargin);
notext2values = ismember('notext2values',varargin);
mergeheaderlines = ismember('mergeheaderlines',varargin);
maketable = ismember('maketable',varargin);

% scan
tab = [];
attributes = struct([]);
for i=1:m % each file
    [pa,na,ex] = fileparts(filename{i});
    if ~isempty(pa) && ~exist(pa,'dir'), error('the directory ''%s'' does not exist.',pa), end
    if (~strcmpi(ex,'.xls')) || (~strcmpi(ex,'.xlsx'))
        if exist([fullfile(pa,na) '.xls'],'file'), ex='.xls'; else ex='.xlsx'; end
    end
    filename{i} = [fullfile(pa,na) ex];
    if ~exist(filename{i},'file'), error('the file ''%s'' does not exist in ''%s''',na,pa), end
    [a,sheets,fmt] = xlsfinfo(filename{i});
    if isempty(a), error('the file ''%s'' in ''%s'' is not recognized',n,p), end
    fileinfo(filename{i})
    if isempty(fmt), disp('No Excel ActiveX server found'), else dispf('Excel format: %s',fmt), end
    if allsheetson, sheetname = sheets; n=length(sheets);
    elseif firstsheeton, sheetname = sheets(1); n=1;
    elseif lastsheeton, sheetname = sheets(end); n=1;
    end
    [localpath,localfile]=fileparts(filename{i}); %#ok<ASGLU>
    filename_clean = replacefromtable(localfile,table);
    if filename_clean(1)<='9', filename_clean = ['file_' filename_clean]; end %#ok<AGROW>
    for j=1:n % each sheet
        sheetname_clean = replacefromtable(sheetname{j},table);
        if sheetname_clean(1)<='9', sheetname_clean = ['sheet_' sheetname_clean]; end %#ok<AGROW>
        if ~isempty(sheetname{j})
            if ~ismember(sheetname{j},sheets), error('the sheet ''%s'' does not exist in ''%s''',sheetname{j},filename{i}); end
            [num,alpha,raw] = xlsread(filename{i},sheetname{j}); %#ok<ASGLU>
        else
            [num,alpha,raw] = xlsread(filename{i}); %#ok<ASGLU>
        end
        if transposeon, raw=raw'; end
        if size(raw,1)>sum(headerlines)
            if length(headerlines)>1
                headers = raw(sum(headerlines(1:end-1))+1:sum(headerlines),:); % header lines
                raw = raw(sum(headerlines)+1:end,:); % raw values
            else
                headers = raw(1:headerlines,:); % header lines
                raw = raw(headerlines+1:end,:); % raw values
                if mergeheaderlines
                    for ih=1:size(headers,2) % added on 19/09/2015
                        headers{1,ih} = regexprep(sprintf(['%s' varlineseparator],headers{:,ih}),varlineseparatorpat,'');
                    end
                    headers = headers(1,:);
                end
            end
        else
            headers = repmat({'unknown'},1,size(raw,2));
        end
        [mraw,nraw] = size(raw);
        if isempty(sheetname{j})
            dispf('[%s][%s] headers:%d table:%dx%d',localfile,sheets{1},headerlines(end),mraw,nraw)
        else
            dispf('[%s][%s] headers:%d table:%dx%d',localfile,sheetname{j},headerlines(end),mraw,nraw)
        end
        if ~iscell(raw), raw = {raw}; end % added OV 04/08/10
        if ~isempty(raw)
            ichar = find(cellfun('isclass',raw,'char')); % index of char cells
            for kw = kwlist, raw(ichar(strcmpi(raw(ichar),kw{1})))={str2double(kw{1})}; end % Conversions 'NaN' 'Inf' '+Inf' '-Inf'
            if ~emptyisnotnan, raw(cellfun('isempty',raw))={'NaN'}; end % replace empty cells by NaNs
            % headers/variables
            var = cell(1,nraw);
            for k=1:nraw
                if ischar(headers{k}), var{k} = replacefromtable(headers{k},table); else var{k} = ''; end
                varlen = length(var{k});
                if varlen>VARMAXLEN, var{k} = var{k}(1:VARMAXLEN); end
                if isempty(var{k}) || nohearderson, var{k} = sprintf('%s%d',varprefix,k); end
                while ismember(var{k},var(1:k-1)), var{k} = [var{k} duplicatesuffix]; end
            end
            [tmp,attr] = deal(cell2struct(repmat({[]},1,nraw),var,2));
            % columns content
            for k=1:nraw
                if all(cellfun('isclass',raw(:,k),'double')) % numerical vector
                    if transposeon, val = [raw{:,k}]; else val = [raw{:,k}]'; end
                else % otherwise cell vector
                    val = raw(:,k);
                    if ~notext2values % automatic conversion of text numbers into values
                        ichar = find(cellfun('isclass',val,'char'));
                        cval = cellfun(@str2double,val(ichar));
                        iconv = ~isnan(cval);
                        cval = cval(iconv);
                        if ~isempty(cval)
                            val(ichar(iconv)) = num2cell(cval); %mat2cell(cval,ones(size(cval)),1);
                            if all(cellfun('isclass',val,'double'))
                                if transposeon, val = [val{:}]; else val = [val{:}]'; end
                            end
                        end
                    end
                end
                tmp.(var{k}) = val;
                attr.(var{k}) = struct('filename',filename,'table',sheetname{j},'title',var{k},'description','');
                attr.(var{k}).description = headers{k}; % prevent the expansion of headers when a cell is used in a struct
                if isnumeric(attr.(var{k}).description), attr.(var{k}).description = ''; end
            end % each column
            %
            if maketable % works only on recent versions of Matlab (2012)
                tmp = struct2table(tmp);
                attrtmp = struct2cell(attr); attrtmp = [attrtmp{:}];
                tmp.Properties.UserData = attrtmp;
                tmp.Properties.Description = sprintf('%s::%s:%s',attrtmp(1).filename,attrtmp(1).table,attrtmp(1).title);
                tmp.Properties.VariableDescriptions = {attrtmp.description};
            end
            % final form
            if m>1 && n>1,      tab.(filename_clean).(sheetname_clean) = tmp;
                                attributes(1).(filename_clean).(sheetname_clean) = attr;
            elseif m==1 && n>1, tab.(sheetname_clean) = tmp;
                                attributes(1).(sheetname_clean) = attr;
            elseif m>1 && n==1, tab.(filename_clean) = tmp;
                                attributes(1).(filename_clean) = attr;
            else                tab = tmp;
                                attributes = attr;
            end
        else
            dispf('\t\t empty sheet (discarded)')
        end
    end % each sheet
end % each file

if nargout>1, attrout = attributes; end