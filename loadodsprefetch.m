function [data,isupdated,attrout]=loadodsprefetch(filename,varargin)
%LOADODSPREFETCH is a common surrogate for LOADODS and XLSTBLREAD enabling and managing prefetch files
% SYNTAX (it uses the same syntax as LOADODS, the syntax is made back compatible to XLSTBLREAD)
%     data = loadodsprefetch(filename [,property1,value1,property2,value2,...])
%     [data,isupdated]=loadodsprefetch(...), isupdated=true if the data has been updated
%
%       filename with .ODS as extension forces the use of LOADODS, filename with .XLS/.XLSX forces the use of XLSTBLREAD
%
%   LIST OF PAIR PROPERTY/VALUE (DEFAULT VALUE)
%         'sheetname' : ''    ... sheet to load (empty for the first sheet, 'all' for all sheets)
%             'blank': NaN    ... value for blank cells
%         'transpose': false  ... true to transpose table/headers (available only with LOADODS)
%            'concat': true   ... true to concat column content (available only with LOADODS)
%           'headers': 1      ... number of header lines
%            'struct': true   ... to transform the table as a structure with headers (available only with LOADODS)
%       'structarray': false  ... to transform the original structure of arrays into a structure array (available only with LOADODS)
%      'forceboolean': false  ... convert to a boolean any column with uniquely 0 and 1 (available only with LOADODS)
%         'maketable': false  ... convert generated structures into tables with attributes (works only on recent versions of Matlab, since 2012)
%   
%    SPECIFIC PROPERTIES TO MANAGE PREFETCH FILES
%           prefetchprefix: prefetch extension (default = 'PREFETCH_')
%             prefetchpath: path of the prefetch file (default=tempdir)            
%               noprefetch: flag to forcce the prefetch to be updated (default=false)
%            prefetchsheet: flag to force a different prefetch for each worksheet (default=false, forced to true if realname is used)
%                 realname: name to be used in the final structure instead of the sheetname 
%                           TIP 1: any empty value will be replaced by the sheet name
%                           TIP 2: if realname is too short, corresponding values of sheetname are used
%
% SEE ALSO: LOADODS, XLSTBLREAD, BYKEYWORDS
%
%
% TWO ADVANCED EXAMPLES FROM THE PROJECT FINDUS
%
% 1) READ MULTIPLE SHEETS (WITH DIFFERENT PREFETCH FILES) AS A TABLE AND RENAME THEM
%    NOTE THAT REALNAME DOES NOT NEEDD TO BE ON THE SAME SIZE (AS THE SUBSITUTION CAN BE PARTIAL
%   db= loadodsprefetch(fullfile(local,datafile),...
%    'sheetname',{'variables','boxes','samples','RHsalts','DVS','TGA','DSC','oven','awmeter','layout'},...
%    'realname', {'vars'     ,''     ,''       ,'RHtable'},...
%    'maketable',true);
%
% 2) READ TEMPERATURES STORED IN DIFFERENT EXCEL FILES (WITH LABVIEW)
%   f = explore('*.xlsx',local,[],'abbreviate'); nf = length(f);
%   for i=1:nf
%       [f(i).data,~,f(i).attr] = loadodsprefetch(fullfile(f(i).path,f(i).file));
%   end

% MS 2.1 - 20/01/12 - INRA\Olivier Vitrac rev. 08/11/15

% Revision history
% 24/01/12 add a comparison based on requested sheetnames
% 26/01/12 use 'case' with argcheck to keep case in folder and file names.
% 26/01/12 fix non-empty loadodsoptions (as a list instead of a structure)
% 12/06/12 fix error message file missing
% 20/09/13 add isupdated
% 08/12/13 force prefetchupdate = true when new spreadsheets are requested
% 17/01/14 force columnwise loadodsoptions
% 07/11/15 manage different prefetch and names for each sheet, first implementation of XLSTBLREAD
% 08/11/15 updated help
% 15/03/16 add xlsm format
%          add 'headerrowindex' property for excel file reading

% default
default = struct(...
    'prefetchprefix','PREFETCH_',...
    'prefetchpath',tempdir,...
    'noprefetch',false,...
    'prefetchsheet',false,...
    'sheetname',[],...
    'realname',[] ...
    );
validchars = '[^a-zA-Z0-9]'; % accepted characters for fields
propertiesODS2XLS = struct(...
          'blank',NaN   ,... value for blank cells
      'transpose',false ,... true to transpose table/headers
        'headers',1     ,... number of header lines
 'headerrowindex',1     ,... index of header line (see xlstblread)
      'noheaders',false,...
  'emptyisnotnan',false,...
  'notext2values',false,...
'mergeheaderlines',false,...
       'maketable',false ...
);

% arg check
if nargin<1, error('one argument is at least required'), end
[options,loadodsoptions] = argcheck(varargin,default,'','case');
if ~isempty(options.realname), options.prefetchsheet = true; end

% manage different prefetch and names for each sheet (added 07/11/2015)
if isempty(options.sheetname)
    options.sheetname='';
elseif options.prefetchsheet
    if ~iscell(options.sheetname), options.sheetname = {options.sheetname}; end
    if isempty(options.realname), options.realname = options.sheetname; end
    if ~iscell(options.realname), options.realname = {options.realname}; end
    nsheets = length(options.sheetname); nreals = length(options.realname);
    options.realname = argpad(options.realname,nsheets,'periodic',options.sheetname(nreals+1:end));
    if nsheets==0 || ismember('all',options.sheetname)
        error('sheetname='''' or sheetname=''all'' are not authorized with the flag prefetchsheet, list all sheets explicitly')
    end
    updtd = false; [data,attributes] =deal(struct([]));
    for i=1:nsheets
        wsk = regexprep(options.sheetname{i},validchars,'');
        wskreal = regexprep(options.realname{i},validchars,'');
        if isempty(wskreal), wskreal = wsk; end
        [data(1).(wskreal),isupdatedtmp,attributes(1).(wskreal)] = ...
            loadodsprefetch(...
            filename,'sheetname',wsk,'realname','','prefetchsheet',false,...
            'prefetchpath',options.prefetchpath,...
            'prefetchprefix',sprintf('%s_%s_',options.prefetchprefix,wsk),...
            loadodsoptions{:} ... all options except sheetname
            );
        updtd = isupdatedtmp || updtd;
    end
    if nsheets==1, data = data.(wskreal); end
    if nargout>1, isupdated = updtd; end
    if nargout>2, attrout = attributes; end
    return    
end

% arg check (cont'ed)
if ~isempty(options.sheetname), loadodsoptions(end+1:end+2) = {'sheetname';options.sheetname}; end % %propagate sheetname
loadodsoptions = loadodsoptions(:);
[~,prefetchfile] = fileparts(filename);
if ~exist(filename,'file'), error('the supplied file ''%s'' does not exist',filename); end
[~,~,ext] = fileparts(filename); ext = lower(ext);
isods = strcmp(ext,'.ods');
if ~isods && ~strcmp(ext,'.xls') && ~strcmp(ext,'.xlsx') && ~strcmp(ext,'.xlsm'); error('only ODS, XLS, XLSX and XLSM files are managed'), end
if isunix && ~isods, error('XLSTBLREAD and XLSREAD requires a WINDOWS environment to run'), end
prefetchfile = fullfile(options.prefetchpath,[options.prefetchprefix prefetchfile '.mat']);
prefetchupdate = false;

% check that prefetch is up to date or match current needs
useprefetch = false; % default behavior
if options.noprefetch
    dispf('LOADODSPREFETCH: noprefetch option is used. The prefetch file is not used')
    prefetchupdate = true;
elseif ~exist(prefetchfile,'file')
    dispf('LOADODSPREFETCH: no prefetchfile detected')
    prefetchupdate = true;
else
    ref = dir(filename);
    pre = dir(prefetchfile);
    load(prefetchfile,'nfo');
    if (ref.datenum<pre.datenum) && (ref.datenum==nfo.datenum) && (ref.bytes==nfo.bytes) %#ok<NODEF> % prefetch up-to-date and same size
        load(prefetchfile,'sheetname');
        if ischar(sheetname) %#ok<NODEF>
            if strcmp(sheetname,'all')
                tmp=load(prefetchfile,'data');
                if isstruct(tmp.data), 
                    sheetname = fieldnames(tmp.data);
                    if ischar(options.sheetname) && strcmp(options.sheetname,'all'), options.sheetname = sheetname; end
                else
                    sheetname = {sheetname};
                end
            else
                sheetname = {sheetname};
            end
        end 
        if ischar(options.sheetname), options.sheetname = {options.sheetname}; end
        if isempty(setdiff(options.sheetname,sheetname))
            useprefetch = true;
        else
            dispf('LOADODSPREFETCH: prefetchfile is not used as additional worksheets are asked to be loaded.')
            prefetchupdate = true;
        end
    else
        dispf('LOADODSPREFETCH: prefetchfile is obsolete')
        prefetchupdate = true;
    end
end

% load data
if useprefetch
    dispf('LOADODSPREFETCH: use the prefetchfile below')
    fileinfo(prefetchfile)
    load(prefetchfile,'data','attributes') % load data
    if length(sheetname)>1 || strcmp(sheetname,'all') 
        data = rmfield(data,setdiff(sheetname,options.sheetname)); %#ok<NODEF> % remove unwanted worksheets
        attributes = rmfield(attributes,setdiff(sheetname,options.sheetname)); %#ok<NODEF>
    end
    fdata = fieldnames(data);
    if length(fdata)==1 && isstruct(data.(fdata{1}))
        data = data.(fdata{1}); attributes = attributes.(fdata{1});
    end
else
    if isods % ODS file
        loadodsoptions = [{filename};loadodsoptions];
        [data,attributes] = loadods(loadodsoptions{:});
    else % XLS/XSLX file
        if ~isempty(loadodsoptions)
            xlstblreadoptions = argcheck(loadodsoptions,[]);
            if isempty(xlstblreadoptions), error('please, replace keywords by pair protery/boolean value'), end
        else
            xlstblreadoptions = loadodsoptions;
        end
        xlstblreadoptions = argcheck(xlstblreadoptions,propertiesODS2XLS,'','keep'); 
        xlstblkeywords = struct(...
            'transpose',xlstblreadoptions.transpose,...
            'noheaders',xlstblreadoptions.noheaders,...
            'emptyisnotnan',~isnan(xlstblreadoptions.blank),...
            'notext2values',xlstblreadoptions.notext2values,...
            'mergeheaderlines',xlstblreadoptions.mergeheaderlines,...
            'maketable',xlstblreadoptions.maketable);
        kw = fieldnames(xlstblkeywords); kw = kw(structfun(@any,xlstblkeywords));
        [data,attributes] = xlstblread(filename,options.sheetname,xlstblreadoptions.headers,'headerrowindex',xlstblreadoptions.headerrowindex,kw{:});
    end
end

% update prefetch file if needed
if prefetchupdate
    nfo = dir(filename); %#ok<NASGU>
    sheetname = options.sheetname; %#ok<NASGU>
    save(prefetchfile,'data','nfo','sheetname','attributes')
    dispf('LOADODSPREFETCH: the prefetch file below has been updated')
    fileinfo(prefetchfile)
end

if nargout>1, isupdated=prefetchupdate; end
if nargout>2, attrout = attributes; end