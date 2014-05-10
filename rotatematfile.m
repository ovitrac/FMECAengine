function rotatematfile(f,freq,n,nofirstpropagate)
%ROTATEMATFILE as rotatelog does but on valid MAT files (corrupted MAT files are not saved)
%   syntax: rotatematfile(f [,freq,n,nofirstpropagate])
%   f = filename
%   freq = frequency in days (default =2)
%   n = number of rotations (default = 4)
%   nofirstpropagate (default=false), use true to save disk space

% ADMIN www v. 2.0 - 12/10/10 - INRA\Olivier Vitrac - rev. 29/12/12


% revision history
% 23/12/12 fix n and help
% 29/12/12 add nofirstpropagate

% default
freq_default = 2;
n_default = 4;
nofirstpropagate_default = false;
null = struct('localname',localname,'date',datestr(now)); %#ok<NASGU>


% arg check
if nargin<1, error('one argument is required'), end
if ~ischar(f), error('f must be a string'), end
[pathstr,fn,ext] = fileparts(f);
if ~strcmpi(ext,'.mat'), error('ONLY MAT files are accepted'), end
if isempty(pathstr), pathstr = pwd; end
if ~exist(f,'file'), error('the file MAT file ''%s%s'' does not exist in ''%s''',fn,ext,pathstr), end
if nargin<2, freq = []; end % every 2 days
if nargin<3, n = []; end % four rotation
if nargin<4, nofirstpropagate = []; end
if isempty(freq), freq = freq_default; end
if isempty(n), n = n_default; end
if isempty(nofirstpropagate), nofirstpropagate = nofirstpropagate_default; end

% build file list
flist = [{f} arrayfun(@(i)regexprep(f,'(.mat)$',sprintf('.%d$1',i)),1:n,'UniformOutput',false)];


% first propagation (if required)
% 1: most recent
% 2: less recent
%..
%n+1 oldest
for i=2:n+1
    if ~exist(flist{i},'file')
        if nofirstpropagate
            save(flist{i},'null')
        else
            back(flist{i-1},flist{i})
        end
    end
end

% retropropagation if files are aged
% snfo = dir(flist{1}); sage = datenum(snfo.date);
% dnfo = dir(flist{2}); dage = datenum(dnfo.date);
% if sage-dage>=freq
if age(flist{1})-age(flist{2})>=freq
    for i=n+1:-1:2
        back(flist{i-1},flist{i});
    end
end

function back(s,d)
%SIMBLE BACKUP function
try %#ok<TRYNC>
    if exist(s,'file')
        data = load(s); %#ok<NASGU>
        save(d,'-struct','data')
    end
end

function a=age(f)
% datenum of file (robust for old Matlab versions)
if exist(f,'file')
    a = dir(f);
    if isfield(a,'datenum'), a=a.datenum; else a=datenum(a.date); end
else
    a = now;
end
