function kin = fmecamerge(fmecadbfile)
%MERGEkinetics merge a fmecadatabase into common kinetics (one per path)
%   kin = fmecamerge(fmecadbfile)
%   [kin,paths] = fmecamerge(...)

% Migration 2.1 (Fmecaengine v0.50) - 08/05/2014 - INRA\Olivier Vitrac - rev.


if ~nargin || ~ischar(fmecadbfile), error('one argument is required, it must be a valid filename'), end
currentpath = rootdir(fmecadbfile);
currentdb = lastdir(fmecadbfile);
if ~exist(fmecadbfile,'file'), error('''%s'' does not exist in ''%s''',currentdb,currentpath), end

% load the databae
dispf('\nFMECAMERGE resuses the following database'), fileinfo(fmecadbfile)
db = load(fmecadbfile);

% build all paths
paths = buildmarkov(db.fmecadb);
npaths = length(paths);

% read all paths
kin = repmat(struct('dbfile','','path','','subsdbfiles','','t',[],'CF',[]),npaths,1);
for ipath = 1:npaths
    steps = paths{ipath}; nsteps = length(steps);
    % read all steps
    currenttime = 0; [CF,tCF] = deal([]); screen = ''; subsdbfiles = cell(nsteps,1); t0 = clock;
    for currentstep = 1:nsteps
        screen = dispb(screen,'\t[path %d/%d] concatenating ''%s'' step %d of %d',ipath,npaths,currentdb,currentstep,nsteps);
        subsdbfiles{currentstep} = fullfile(currentpath,steps{currentstep});
        currentres = load(subsdbfiles{currentstep}); % load r
        if currentstep==1
            validtimes = ((currenttime+currentres.r.t*currentres.r.timebase)>=currenttime) & ...
                         (currentres.r.t*currentres.r.timebase<=db.fmecadb.(steps{currentstep}).t);
        else
            validtimes = ((currenttime+currentres.r.t*currentres.r.timebase)>currenttime) & ...
                         (currentres.r.t*currentres.r.timebase<=db.fmecadb.(steps{currentstep}).t);
        end
        tCF = [tCF;currenttime+currentres.r.t(validtimes)*currentres.r.timebase]; %#ok<AGROW>
        CF = [CF;currentres.r.CF(validtimes)]; %#ok<AGROW>
        currenttime = tCF(end);
    end % next current step
    kin(ipath) = struct('dbfile',fmecadbfile,'path',{steps'},'subsdbfiles',{subsdbfiles},'t',tCF,'CF',CF);
    dispb(screen,'\t[path %d/%d] completed %d steps in %0.4g s',ipath, npaths,nsteps,etime(clock,t0));
end % next path