function fmecacmd(fmecamainfile,fmecaoptionfile)
% Syntax:
%   fmecacmd myfmeca.ods
%   fmecacmd myfmeca.ods setup.xml
% INRA - Mai Nguyen, Olivier Vitrac
% history - 18/06/2018

% compilation directive 
% mcc -mv fmecacmd.m -a Dpiringer.m

defaultfmeca = struct(...
    'nmesh', 200,...   % refined 500
    'reltol', 1e-4,... % refined 1e-6
    'abstol', 1e-4,... % refined 1e-6
    'initialstep',1e-8,...
    'maxstep',1e3,...
    'maxorder',2,...
    'inputpath','inputs',...
    'outputpath','outputs',...
    'fmecasheetname','sim'...
     );
 
% arg check
if nargin < 1, error('Please supply an input file'), end
if nargin < 2, fmecaoptionfile = ''; end

% check fmecamainfile
if ~ischar(fmecamainfile), error('the first parameter should be a string'), end
[fpath, fmecaname, fext] = fileparts(fmecamainfile);
if ~strcmpi(fext,'.ods'), error('The input file should be an ODS file'), end
if isempty(fpath), fpath = pwd; end
if ~exist(fpath,'dir')
    error('The supplied folder ''%s'' of input file does not exist',fpath)
end
if ~exist(fullfile(fpath,[fmecaname fext]),'file')
    error('The supplied input file ''%s'' does not exist in ''%s''',[fmecaname fext],fpath)
end
local = fpath; 
fmecamainfile_original = fullfile(fpath,[fmecaname fext]);

% check fmecaootionfile
if ~isempty(fmecaoptionfile)
    if ~ischar(fmecaoptionfile), error('the second parameter should be a string'), end
    [fpath, fname, fext] = fileparts(fmecaoptionfile);
    if isempty(fpath)
        fpath = pwd;
        if ~exist(fullfile(fpath,[fname fext]),'file')
            fpath = local;
        end
    end
    if ~strcmpi(fext,'.xml'), error('The option file should be an XML file'), end
    if ~exist(fpath,'dir')
        error('The supplied folder ''%s'' of option file does not exist',fpath)
    end
    if ~exist(fullfile(fpath,[fname fext]),'file')
        error('The supplied option file ''%s'' does not exist in ''%s''',[fname fext],fpath)
    end
    options = loadSFPP3(fmecaoptionfile);
    if isfield(options,'options'), error('OPTIONS fied does not exist in ''%s''',fmecaoptionfile), end
    options = argcheck(options.options,defaultfmeca);
    saveSFPP3(fullfile(local,fname),{options},{'options'})
else
    options = defaultfmeca;
    fmecaoptionfile = fullfile(local,'options');
    saveSFPP3(fmecaoptionfile,{options},{'options'})
end

% duplicate input file in %local%\inputs
fmecamainfile = fullfile(local,options.inputpath,[fmecaname '.ods']);
if ~exist(fullfile(local,options.inputpath),'dir')
    [sts,msg] = mkdir(fullfile(local,options.inputpath));
    if ~sts, error(msg), end
end
[sts, msg] = copyfile(fmecamainfile_original,fmecamainfile);
if ~sts, error(msg), end

% reference simulations
solveropt = odeset('RelTol',options.reltol,'AbsTol',options.abstol,'initialstep',options.initialstep,...
                 'Maxstep',options.maxstep,'Maxorder',options.maxorder);
fmecadef = struct( ...
    'local',local,...
    'fmecasheetname',options.fmecasheetname,... % datasheet to be read of ods file
    'fmecamainfile',[fmecaname '.ods'],... FMECAengine assumes always local\inputs
    'database',struct([]),... 
    'inputpath',options.inputpath,...
    'outputpath',options.outputpath, ... 
    'fmecadbfile','fmecamodal.mat',...
    'nmesh',options.nmesh,...
    'options',solveropt,...
    'nograph','true'...
    );
fmecaengine(fmecadef)