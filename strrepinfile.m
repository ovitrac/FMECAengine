function ret=strrepinfile(infile,outfile,rule,varargin)
%STRREPINFILE Replace string with another in file (regular expression are accepted)
%   syntax1: strrepinfile(infile,outfile,rule [,'property',value])
%             infile: input filename
%            outfile: output filename (if empty, outfile=infile)
%               rule: set of simple replacement rules coded within a structure
%                   rule.texttoreplace = replacingtext
%             recognized properties include
%               'suffix' : suffix string to append to rule.texttoreplace (default = '')
%               'prefix': prefix string to prepend to rule.texttoreplace (default = '')
%               'backupsuffix': suffix added to the backup of outfile if it exist (default = '~')
%
%   syntax2: strrepinfile(infile,outfile,stringstofind,replacingstrings [,'backupsuffix',string])
%             infile: input filename
%            outfile: output filename (if empty, outfile=infile)
%      stringstofind: string or cell array of strings coding for regular expression
%   replacingstrings: string or cell array of strings coding for regular expression
%
%   option: ret = strrepinfile(...) returns 0 if successful or -1 if not
%
%   Note: syntax1 uses strrep
%         syntax2 uses regexprep, protect characters if required (with regexptranslate)
%
%
%   Example: replace all occurrences of $$lookfor1$$ and $$lookfor2$$ in an HTML file 
%            by replacementstring1 and replacementstring2 respectively
%       local = 'C:\Data\Olivier\Audrey_Goujon\www\SFPD\questionnaire'
%       infile = 'source.html';
%       outfile = 'destination.html';
%       clear rules
%       rules.lookfor1 = 'replacementstring1';
%       rules.lookfor2 = 'replacementstring2';
%       strrepinfile(fullfile(local,infile),fullfile(local,outfile),rules,'prefix','$$','suffix','$$');

% Migration 2.0 - 15/04/2011  INRA\Olivier Vitrac - rev.

% Revision history

% default
options_default = struct(...
    'suffix','',...
    'prefix','',...
    'backupsuffix','~' ...
    );
if isunix, eol = '\n'; else eol = '\r\n'; end

% argcheck
if nargin<3, error('three arguments are required'), end
if ~ischar(infile), error('the first argument must be a filename'), end
if ~isempty(outfile) && ~ischar(outfile), error('the second argument must be a filename'), end
if isempty(outfile), outfile = infile; end
if ~exist(infile,'file')
    [filestr,pathstr] = lastdir(infile);
    if isempty(infile), pathstr = pwd; end
    error('the file ''%s'' does not exist in ''%s''',filestr,pathstr)
end
if isstruct(rule)
    syntax1 = true;
    options = argcheck(varargin,options_default);
else
    if nargin<4, error('four arguments are required'); end
    syntax1 = false;
    stringstofind = rule;
    replacingstrings = varargin{1};
    if ~iscell(stringstofind), stringstofind = {stringstofind}; end
    if ~iscell(replacingstrings), replacingstrings = {replacingstrings}; end
    if ~iscellstr(stringstofind), error('stringstofind must be a string or a cell array of strings'); end
    if ~iscellstr(replacingstrings), error('replacingstrings must be a string or a cell array of strings'); end
    nstringstofind = length(stringstofind);
    nreplacingstrings = length(replacingstrings);
    if (nstringstofind>1) && (nreplacingstrings>1) && (nstringstofind~=nreplacingstrings)
        error('the sizes of stringstofind and replacingstrings are not compatible')
    end
    options = argcheck(varargin(2:end),options_default);
end

% load file
fid = fopen(infile,'r');
txt = uncell(textscan(fid,'%s','delimiter','\n'));
fclose(fid);


% replacements
if syntax1
    for f=fieldnames(rule)'
        txt = strrep(txt,sprintf('%s%s%s',options.prefix,f{1},options.suffix),rule(1).(f{1}));
    end
else
    txt = regexprep(txt,stringstofind,replacingstrings);
end

% backup if required
if exist(outfile,'file')
    success = movefile(outfile,sprintf('%s%s',outfile,options.backupsuffix));
    if ~success, error('unable to backup the existing output ''%s''',outfile); end
end

% save
fid = fopen(outfile,'w');
fprintf(fid,['%s' eol],txt{:});
tmp = fclose(fid);
if tmp, error('unable to write ''%s''',outfile); end

% output
if nargout, ret = tmp; end