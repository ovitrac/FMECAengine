function data=loadodsprefetch(filename,varargin)
%LOADODSPREFETCH loadods surrogate to use/manage prefetch files when they exist
%      data = loadodsprefetch(...)
%         It uses the same syntax as loadods.
%         Additional properties/values are:
%           prefetchprefix: prefetch extension (default = 'PREFETCH_')
%             prefetchpath: path of the prefetch file (default=tempdir)            
%               noprefetch: flag to forcce the prefetch to be updated (default=false);


% MS 2.1 - 20/01/12 - INRA\Olivier Vitrac rev.

% default
default = struct(...
    'prefetchprefix','PREFETCH_',...
    'prefetchpath',tempdir,...
    'noprefetch',false ...
    );

% arg check
if nargin<1, error('one argument is at least required'), end
[options,loadodsoptions] = argcheck(varargin,default);
[~,prefetchfile] = fileparts(filename);
if ~exist(filename,'file'), error('the supplied file ''%s'' does not exist'); end
prefetchfile = fullfile(options.prefetchpath,[options.prefetchprefix prefetchfile '.mat']);
prefetchupdate = false;

% load section
if options.noprefetch
    dispf('LOADODSPREFETCH: noprefetch option is used. The prefetch file is not used')
    data = loadods(filename,loadodsoptions);
    prefetchupdate = true;
elseif ~exist(prefetchfile,'file')
    dispf('LOADODSPREFETCH: no prefetchfile detected')
    data = loadods(filename,loadodsoptions);
    prefetchupdate = true;
else
    ref = dir(filename);
    pre = dir(prefetchfile);
    load(prefetchfile,'nfo');
    if (ref.datenum<pre.datenum) && (ref.datenum==nfo.datenum) && (ref.bytes==nfo.bytes) %#ok<NODEF> % prefetch up-to-date and same size
        dispf('LOADODSPREFETCH: use the prefetchfile below')
        fileinfo(prefetchfile)
        load(prefetchfile,'data') % load data
    else
        dispf('LOADODSPREFETCH: prefetchfile is obsolete')
        data = loadods(filename,loadodsoptions);
        prefetchupdate = true;
    end
end

% update prefetch file if needed
if prefetchupdate
    nfo = dir(filename);
    save(prefetchfile,'data','nfo')
    dispf('LOADODSPREFETCH: the prefetch file below has been updated')
    fileinfo(prefetchfile)
end
