function db=loadfmecaenginedb(varargin)
%LOADFMECAENGINEDB load an FMECAengine database (e.g. database of substances) using a prefetch file (if needed)
%   syntax: db = loadfmecaenginedb;
%           db = loadfmecaenginedb();
%           db = loadfmecaenginedb(property1,value1,...,keyword);
%   Property/value(optional)
%       dbfile: filename to an ODS database (default = 'substancedb.ods')
%       dbpath: full path (default = fullfile(find_path_toolbox('migration'),'database'))
%    sheetname: to be used by loadods (default='all'), sheets to load
%      headers: to be used by loadods (default=1), number of header rows
% prefetchfile: filename of the prefetch file (assumed to be locatedd in dbpath, default = PREFETCH_dbfile.mat)
%   Keyword (optional)
%       noprefetch: flag to avoid the use of prefetch (the prefetch will be also regenerated)
%       NB: noprefetch can be also used as a property as 'nopprefetch',noprefetchflag
%
%   Example: db = loadfmecaenginedb('dbfile','mydatabase','sheetname',{'substance' 'polymer'})
%
%   Note: since LOADODSPREFETCH, this function tends to be depreciated.
%         To keep retrocompatibility LAODODS has been internally replaced by LOADODSPREFETCH
%
%   See also: loadods, loadodsprefetch, fmecaengine

% MIGRATION 2.0 - 20/07/11 - INRA\Olivier Vitrac - rev. 08/02/12

% Revision history
% 24/08/11 fix help, enables noprefetch as a property
% 25/08/11 fix case
% 08/02/12 loadodsprefetch is used internally instead of loadods

% default
prop_default = struct( ...
    'dbfile','substancedb.ods',...
    'dbpath',fullfile(find_path_toolbox('migration'),'database'),...
    'prefetchfile','',...
    'sheetname','all',...
    'headers',1,...
    'noprefetch',0 ...
    );
kwlist = 'noprefetch';

% argcheck
prop = argcheck(varargin,prop_default,kwlist,'property','case');
if isempty(prop.prefetchfile)
    prop.prefetchfile = sprintf('PREFETCH_%s',regexprep(prop.dbfile,'.ods$','.mat','ignorecase'));
end

% load prefetch or ods data
if ~prop.noprefetch && exist(fullfile(prop.dbpath,prop.prefetchfile),'file')
    dispf('Load prefetchfile...')
    fileinfo(fullfile(prop.dbpath,prop.prefetchfile))
    load(fullfile(prop.dbpath,prop.prefetchfile))
else
    db = loadodsprefetch(fullfile(prop.dbpath,prop.dbfile),'sheetname',prop.sheetname,'headers',prop.headers);
    save(fullfile(prop.dbpath,prop.prefetchfile),'db')
    dispf('...prefetch file updated')
    fileinfo(fullfile(prop.dbpath,prop.prefetchfile))
end