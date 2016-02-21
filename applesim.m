function Rout = applesim(varargin)
%APPLESIM simulates the migration in and through apple (idealized as a sphere)
%   Default object constructor:
%       S = applesim();
%       S = applesim('initialize','property1',value1,'property2',value2,...)
%   Simulation launcher using the input definition S
%       R = applesim(S)
%   Simulation launcher using the input definition S along with user overrides
%   (note when a property is defined twice or more, only the first definition is kept)
%       R = applesim('property1',value1,'property2',value2,...,S)
%   Advanced use:
%       R = applesim('keyword1','keyword2',...,'property1',value1,'property2',value2,...,S)
%
%  List of pair property/value (property name: default value [meaning, unit])
%  NOTE THAT OPTIONS AND PLOTOPTIONS ARE STRUCTURE (use struct('property1',value1,'property2',value2,...) to set them
%  NOTE THAT EMPTY VALUES (for scalar and structure properties) WILL FORCE INHERITANCE FROM DEFAULT VALUES
%
%              t: [0 logspace(-3,log10(4*24*3600),1000)] [time in seconds]
%          lskin: 1.0000e-03 [skin thickness in m]
%         lflesh: 0.0390     [flesh thickness in m]
%          kskin: 4          [Henry like coefficient for the skin in -]
%         kflesh: 1          [Henry like coefficient for the flesh in -]
%          Dskin: 1.0000e-12 [Diffusion coefficient for the skin in m2/s]
%         Dflesh: 1.0000e-11 [Diffusion coefficient for the flesh in m2/s]
%           Cinf: 1.0000e-04 [external concentration in -]
%     Vinf2apple: 10         [volume external liquid-to-apple ratio in -]
%          nmesh: 200        [number of finite volumes]
%       nmeshmin: 50         [minimum number of finite volumes in each compartment]
%             Bi: 1000       [mass Biot number in -]
%    constructor: @senspatankarCgeometry [engine constructor]
%         engine: @(F) senspatankarCgeometry(senspatankar_wrapper(F)) [simulation engine]
%        options: [1x1 struct] [see odeset]
%    plotoptions: structure setting printing options with fields
%               figname: sprintf('%s_%s',mfilename,datestr(now,'YYYY-MM-DD_HH-mm-ss'))  [figure name, to be used for printing]
%             nprofiles: 12                          [number of concentration profiles to be displayed]
% nprofilesmaxforlegend: 20                          [maximim number of concentration profiles used for legending]
%              colormap: @jet
%                 color: structure with fields (missing colors will be assigned with default ones)
%                              skin: 'DarkOrange'
%                             flesh: 'Chartreuse'
%                             apple: 'DarkMagenta'        
%         paperposition: [0.6345 0.1774 20.3046 29]
%              position: [107 46 845 950]
%                 lunit: 'mm'                        [any length unit recognized by FMECAunit]
%                 tunit: 'day'                       [any time unit recognized by FMECAunit]
%
%  TIPS
%   For fitting, the property 'eval' enables an evaluation of particular fields of the output (see FITTING CASE STUDIES)
%           eval: @(R)... anonymous function depending on R, where R is the output of applesim
%                empty value forces R to be returned (equivalent to @(R) R)
%   For simulating samples of skin or of flesh, set any of the property D, k or l to NaN
%           Example: setting lflesh to NaN forces the simulation of a piece of skin in Cartesian coordinates
%                    setting lskin to NaN forces the simulation of a piece of flesh in Cartesian coordinates
%
%   Keywords:
%       'initialize' returns the default inputs   [example: S = applesim('initialize') ]
%              'run' launches simulation          [example: R = applesim('run') ]
%             'plot' plots the profiles and kinetics
%      'noverbosity' remove verbose messages
%
%  OUTPUT: R = structure with fields (with eval empty)
%        nfo: traceability structure with fields (engine, version, host, user, date)
%     inputs: copy of S used as inputs
%        sim: inputs as sent to the simulation engine
%        raw: raw results returned by the simulation engine
%       skin: skin results
%      flesh: flesh results
%      apple: apple results
%  OPUTPUT with eval active
%       output = S.eval(R)
%
%  Presentation of results for skin, flesh and apple
%        t: nt x 1 vector of real times [unit: s]
%        x: nx x 1 vector of positions counted from the external surface (depths) [m]
%        r: nx x 1 vector of radial positions  [m]
%        C: nx x nt array of concentration profiles [same units as Cinf]
%     Crav: nx x nt array of cumulated concentration profiles  [same units as Cinf]
%      Cav: 1 x nt  [same units as Cinf]
%
%   BASIC EXAMPLE (full apple)
%{
       S = applesim();
       R = applesim(S);
%}
%   EXAMPLE FOR PLOTTING
%{
       plotopt = struct('tunit','hour','lunit','cm');
       R = applesim('plotoptions',plotopt,S,'plot');
%}
%   EXAMPLE FOR STUDYING THE MIGRATION TROUGH A SAMPLE OF SKIN WITHOUT THE FLESH
%{
     R = applesim('Dflesh',NaN,S,'plot'); 
%}
%   EXAMPLE FOR STUDYING THE MIGRATION TROUGH A SAMPLE OF FLESH WITHOUT THE SKIN
%{
     S = applesim('initialize','lskin',NaN,...
                  'lflesh',2.5e-3,...
                  'Dflesh',5e-11,...
                  'kflesh',2.5,...
                  't',linspace(0,sqrt(11*24*3600),1000).^2,...
                  'Vinf2apple',1e3);
     R = applesim('Dflesh',4e-11,S,'plot'); 
%}
%   CASE STUDY (1) FITTING one parameter Dskin
%{
       S = applesim('initialize','Dskin',1.5e-12); %value to guess Dskin
       R = applesim(S,'plot');
       tmeasurements = (0:8:3.5*24)'*3600; % every 8 hours for 3.5 days
       noiselevel = 0.05;
       measurements = struct('t',tmeasurements,'Cav',interp1(R.apple.t,R.apple.Cav,tmeasurements,'pchip').*(1+noiselevel/1.96*randn(size(tmeasurements))));
       figure, plot(R.apple.t,R.apple.Cav,'b-',measurements.t,measurements.Cav,'ro')
       % criterion
       Cprediction = @(Dskin) applesim('eval',@(R) interp1(R.apple.t,R.apple.Cav,tmeasurements,'pchip'),'Dskin',Dskin,S);
       %Cprediction(1e-11) -> gives the concentration profile for Dskin=1e-11 m2/s
       crit = @(logDskin) norm(measurements.Cav-Cprediction(exp(logDskin))); % we fit log(Dskin)
       % minimization of the criterion
       optoptions = optimset('display','iter','FunValCheck','on','MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6);
       [logDopt,err] = fminsearch(crit,log(S.Dflesh),optoptions);
       Dopt = exp(logDopt)
       Ropt = applesim('Dskin',Dopt,S);
       hold on, plot(Ropt.apple.t,Ropt.apple.Cav,'r--')
%}
%   CASE STUDY (2) FITTING simultaneously two parameters: Dskin and kskin
%{
       S = applesim('initialize','Dskin',1.5e-12,'ksin',4.8); %values to guess Dskin and kskin
       R = applesim(S);
       tmeasurements = (0:8:3.5*24)'*3600; % every 8 hours for 3.5 days
       noiselevel = 0.05;
       measurements = struct('t',tmeasurements,'Cav',interp1(R.apple.t,R.apple.Cav,tmeasurements,'pchip').*(1+noiselevel/1.96*randn(size(tmeasurements))));
       figure, plot(R.apple.t,R.apple.Cav,'b-',measurements.t,measurements.Cav,'ro')
       Cprediction = @(Dskin,kskin) applesim('eval',@(R) interp1(R.apple.t,R.apple.Cav,tmeasurements,'pchip'),'Dskin',Dskin,'kskin',kskin,S);
       crit = @(A) norm(measurements.Cav-Cprediction(exp(A(1)),A(2))); % we fit A = [log(Dskin) kskin]
       optoptions = optimset('display','iter','FunValCheck','on','MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6);
       [Aopt,err] = fminsearch(crit,[log(S.Dflesh) S.kflesh],optoptions);
       Dopt = exp(Aopt(1)), kopt = Aopt(2)
       Ropt = applesim('Dskin',Dopt,'kskin',kopt,S);
       hold on, plot(Ropt.apple.t,Ropt.apple.Cav,'r--')

%}
%   CASE STUDY (3) FITTING Dflesh from a sample experiment
%{

          %% get the data
        local = 'C:\FMECAengine-master\data';
        datafile = 'WEAC_partition_coefficient_experiment.xlsx';
        outputfolder = local;
        T = loadodsprefetch(fullfile(local,datafile),'maketable',true,'prefetchpath',local);
        T.data = cleantable(T.data,'none'); % remove inconsistencies in table and replace NaN/empty values in categorical columns by 'none'
        fleshhold = bykeywords(T.data,'type','flesh');
        fleshkeep = fleshhold.kgSrkgapple(4:21);
        tmeasurements = fleshhold.tdaysinSrsolution(4:21).*(3600*24);

          %% format the data (remove surface adhesion (day 0 stuff), remove day 11 outlier
        adhese = zeros(length(fleshkeep),1);
            for i = 1:3;
            adhese(i) = fleshkeep(i);
            end
            for i = 4:length(fleshkeep);
            adhese (i) = mean(fleshkeep(1:3));
            end
        fleshwo = fleshkeep - adhese;
        fleshdata = fleshwo;
        fleshdata(17) = [];
        tw = tmeasurements;
        tw(17) = [];

        S = applesim('initialize','Dflesh',1.1e-13,'kflesh',1.06,'Cinf',0.05,'lflesh',5e-4,'Dskin',NaN); %value to guess Dflesh
        R = applesim(S); 
          %% criterion
       Cprediction = @(Dflesh) applesim('eval',@(R) interp1(R.flesh.t,R.flesh.Cav,tw,'pchip'),'Dflesh',Dflesh,'Dskin',NaN,S);
       crit = @(logDflesh) norm(fleshdata-Cprediction(exp(logDflesh))); % we fit log(Dflesh)

          %% check values using plot beforehand to help initialize to appropriate values
        figure, plot(tw/(24*3600),fleshdata,'o',tw/(24*3600),Cprediction(1.1e-13),'d')
        

          %% minimization of the criterion
       optoptions = optimset('display','iter','FunValCheck','on','MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6);
       [logDopt,err] = fminsearch(crit,log(S.Dflesh),optoptions);
       Dopt = exp(logDopt)
       Ropt = applesim('Dflesh',Dopt,'Dskin',NaN,S);
        figure, plot(tw/(24*3600),fleshdata,'o',tw/(24*3600),Cprediction(Dopt),'d')
       hold on, plot(Ropt.flesh.t/(24*3600),Ropt.flesh.Cav,'r--')

%}
% CASE STUDY (4) FITTING Dflesh and kflesh from a sample experiment
%{
         %% get the data
        local = 'C:\FMECAengine-master\data';
        datafile = 'WEAC_partition_coefficient_experiment.xlsx';
        outputfolder = local;
        T = loadodsprefetch(fullfile(local,datafile),'maketable',true,'prefetchpath',local);
        T.data = cleantable(T.data,'none'); % remove inconsistencies in table and replace NaN/empty values in categorical columns by 'none'
        fleshhold = bykeywords(T.data,'type','flesh');
        fleshkeep = fleshhold.kgSrkgapple(4:21);
        tmeasurements = fleshhold.tdaysinSrsolution(4:21).*(3600*24);

          %% format the data (remove surface adhesion (day 0 stuff), remove day 11 outlier
        adhese = zeros(length(fleshkeep),1);
            for i = 1:3;
            adhese(i) = fleshkeep(i);
            end
            for i = 4:length(fleshkeep);
            adhese (i) = mean(fleshkeep(1:3));
            end
        fleshwo = fleshkeep - adhese;
        fleshdata = fleshwo;
        fleshdata(17) = [];
        tw = tmeasurements;
        tw(17) = [];
        S = applesim('initialize','Dflesh',1.1e-12,'kflesh',0.67,'Cinf',0.05,'lflesh',5e-4,'Dskin',NaN); %value to guess Dflesh and kflesh
        R = applesim(S);
       Cprediction = @(Dflesh,kflesh) applesim('eval',@(R) interp1(R.flesh.t,R.flesh.Cav,tw,'pchip'),'Dflesh',Dflesh,'kflesh',kflesh,'Dskin',NaN,S);
       crit = @(A) norm(fleshdata-Cprediction(exp(A(1)),A(2))); % we fit A = [log(Dskin) kskin]
       optoptions = optimset('display','iter','FunValCheck','on','MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6);
       [Aopt,err] = fminsearch(crit,[log(S.Dflesh) S.kflesh],optoptions);
       Dopt = exp(Aopt(1)), kopt = Aopt(2)
       Ropt = applesim('Dflesh',Dopt,'kflesh',kopt,'Dskin',NaN,S);
       figure, plot(tw/(24*3600),fleshdata,'o',tw/(24*3600),Cprediction(Dopt,kopt),'d')
       hold on, plot(Ropt.flesh.t/(24*3600),Ropt.flesh.Cav,'r--')

%}
%
% See also: senspatankarCgeometry, senspatankar_wrapper
%
%
% RADIOMIG v. 0.31 - 15/10/2015 - INRA\Olivier Vitrac, FDA\Danielle Larese - rev. 20/02/2016


% Revision history
% 15/10/2015 alpha release based on script apple_gen2
% 16/10/2015 add plot capabilities
% 17/10/2015 add plotoptions and its parsing, updated help, release candidate
% 17/10/2015 add fitting examples (version 0.2)
% 19/10/2015 add the ability to simulate a sample of skin or flesh (version 0.3)
% 20/10/2015 add personalized colors and plot capabilities for sample
% 23/10/2015 fix time unit in legend, default time unit is set to 'days' instead of 'day'
% 01/12/2015 force F.C0 = 0; if issample (instead of F.C0 = [0 0])
% 09/12/2015 add 'noverbosity'
% 10/12/2015  add case-study 4
% 19/02/2016 fix nmesh, nmeshmin, options
% 20/02/2016 code cleaning to supress output

%% change the version to match the help of this function
codeversion = 0.306;

%% Default
% list of recognized keywords
keywords = {'initialize' 'run','plot','noverbosity'};
% properties to be used for plotting
defaultcolors = struct('skin','DarkOrange','flesh','Chartreuse','apple','DarkMagenta');
plotoptions = struct('figname',sprintf('%s_%s',mfilename,datestr(now,'YYYY-MM-DD_HH-mm-ss')),...
                     'nprofiles',12,'nprofilesmaxforlegend',20,'colormap',@jet,'color',defaultcolors,...
                     'paperposition',[0.6345    0.1774   20.3046   29],'position',[107    46   845   950],...
                     'lunit','mm','tunit','days');
% ODE15s solver parameters (without 'BDF' on, it is a NDF one's, please no order larger than 2 for BI larger than 10)
options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Initialstep',1e-8,'Maxstep',.01,'Maxorder',2);
% default arguments for simulation, solver, plotting
days = 24*3600; % days in s
default = struct(...
    't',[0 logspace(-3,log10(4*days),1000)],...
    'lskin',1e-3,...
    'lflesh',39e-3,...
    'kskin',4,...
    'kflesh',1,...
    'Dskin',1e-12,...
    'Dflesh',1e-11,...
    'Cinf',100e-6,...
    'Vinf2apple',10,...
    'Bi',1e3,...
    'nmesh',200,...
    'nmeshmin',50,...
    'constructor',@senspatankarCgeometry,...
    'engine',@(F) senspatankarCgeometry(senspatankar_wrapper(F)),...
    'options',options,...
    'plotoptions',plotoptions,...
    'eval',[] ...
        );
% traceability
t0 = clock;
nfo = @(trun) struct('engine',mfilename,'version',codeversion,'host',localname,'user',username,'date',datestr(now),'elapsedtime',etime(trun,t0));
   
%% arg check (please note that the parsing of options and plotoptions is done separately)
if nargin<1, Rout = applesim('initialize'); return, end
[otmp,remain] = argcheck(varargin,struct('options',options,'plotoptions',plotoptions),keywords,'nostructexpand'); % extract options and keywords
otmp.options = argcheck(otmp.options,options); % propagate missing options
otmp.plotoptions = argcheck(otmp.plotoptions,plotoptions); % propagate missing plot options
otmp.plotoptions.color = argcheck(otmp.plotoptions.color,defaultcolors);
o = argcheck(remain,default);   % parse remaining arguments/properties versus default
o = argcheck(otmp,o,'','keep'); % propagate properties

%% check whether a sample is required instead of the whole apple
isskinmissing  =  isnan(o.lskin) || isnan(o.Dskin) || isnan(o.kskin);
isfleshmissing  =  isnan(o.lflesh) || isnan(o.Dflesh) || isnan(o.kflesh);
if isskinmissing  && isfleshmissing
    error('Only skin or flesh can be missing but not all parts')
elseif isfleshmissing 
    issample = true;
    sample = 'skin';    
elseif isskinmissing
    issample = true;
    sample = 'flesh';
else
    issample = false;
    sample = '';
end    
    
%% initialize
if o.initialize
    Rout = o;
    if ~o.run, return, end
end

%% Assembling inputs for senspatankarCgeometry() - other engines can be implemented via engine
F = o.constructor();
F.Bi = o.Bi;
F.t = o.t;
F.CF0 = o.Cinf; % (au) Concentration outside (index 0 = surrounding medium)
F.L = 1/o.Vinf2apple; % Dilution factor (global)
F.C0 = [0 0]; % (au) Initial concentrations (index 1 = skin, index 2 = core)
F.nmesh = o.nmesh;
F.nmeshmin = o.nmeshmin;
F.options = o.options;

if issample
    F.C0 = 0;
    F.ngeometry = 0;
    for prop = {'D' 'k' 'l'}
        F.(prop{1}) = o.(sprintf('%s%s',prop{1},sample));
    end
else
    F.ngeometry = 2;
    F.l = [o.lskin o.lflesh];
    F.D = [o.Dskin o.Dflesh];
    F.k = [o.kskin o.kflesh];
end

%% Launch simulation
R = o.engine(F);

%% Mass balance
t = R.tC*R.timebase;

if issample % Simple mass balance for sample (skin or flesh in cartesian coordinates)

    C = R.Cx';
    x = R.x*R.F.lengthscale;
    Ctot = trapz(x,C)/F.l;
    Cxtot = cumtrapz(x,C) ./ x(:,ones(1,size(C,2)));

else % Mass balances for full apple

    % mass balance for the flesh/parenchyma
    x = R.x*R.F.lengthscale;
    isflesh = x>=F.l(1);
    Cflesh = R.Cx(:,isflesh)';
    xflesh = x(isflesh);
    rflesh = sum(F.l)-xflesh;
    Cflesh_av = trapz(rflesh,(rflesh(:,ones(1,size(Cflesh,2))).^2).*Cflesh) ./ trapz(rflesh,rflesh.^2);
    % mass balance for the skin
    isskin = x<=F.l(1);
    Cskin = R.Cx(:,isskin)';
    xskin = x(isskin);
    rskin = sum(F.l)-xskin;
    Cskin_av = trapz(rskin,(rskin(:,ones(1,size(Cskin,2))).^2).*Cskin) ./ trapz(rskin,rskin.^2);

    % cumulated profile (added INRA\OV - 24/09/2015)
    C = flipud(R.Cx');
    r = flipud(sum(F.l) - R.x*R.F.lengthscale);
    Ctot = trapz(r,(r(:,ones(1,size(C,2))).^2).*C) ./ trapz(r,r.^2);
    Crtot = Ctot(ones(length(r),1),:) - cumtrapz(r,(r(:,ones(1,size(C,2))).^2).*C) ./ cumtrapz(r,r(:,ones(1,size(C,2))).^2);
    
end % issample
    
%% Final result
if issample % simple simulation
    
    Rout = struct(...
        'nfo',nfo(clock),...
        'coordinates','cartesian',...
        'inputs',o,...
        'sim',F,...
        'raw',R,...
        sample,struct(...
            't',t,...
            'x',x,...
            'C',C,...
            'Cxav',Cxtot,...
            'Cav',Ctot )...
        );
    
else % full simulation
    
    Rout = struct(...
        'nfo',nfo(clock),...
        'coordinates','spherical',...
        'inputs',o,...
        'sim',F,...
        'raw',R,...
        'skin',struct(...
            't',t,...
            'is',isskin,...
            'x',xskin,...
            'r',rskin,...
            'C',Cskin,...
            'Cav',Cskin_av ),...
        'flesh',struct(...
            't',t,...
            'is',isflesh,...
            'x',xflesh,...
            'r',rflesh,...
            'C',Cflesh,...
            'Cav',Cflesh_av ),...
        'apple',struct(...
            't',t,...
            'x',x,...
            'r',r,...
            'C',C,...
            'Crav',Crtot,...
            'Cav',Ctot ) ...
            );
        
end
%% eval if required
if ~isempty(o.eval)
    if isa(o.eval,'function_handle')
        Rout = o.eval(Rout);
        return
    else
        error('eval must define an anonymous function')
    end
end
    
%% Do plot if requested
if o.plot
    
    % definitions
    formatfig(figure,'figname',o.plotoptions.figname,'paperposition',o.plotoptions.paperposition,'position',o.plotoptions.position)
    tfit = linspace(0,sqrt(t(end)),o.plotoptions.nprofiles)'.^2; % square root time scale (to get a nice display)
    nleg = min(o.plotoptions.nprofiles,o.plotoptions.nprofilesmaxforlegend);
    col = o.plotoptions.colormap(o.plotoptions.nprofiles); istoolight = mean(col,2)>0.8; col(istoolight) = col(istoolight)-0.2;
    tunit = FMECAunit('t',sprintf('1 %s',o.plotoptions.tunit)); tunitleg = o.plotoptions.tunit;
    lunit = FMECAunit('l',sprintf('1 %s',o.plotoptions.lunit));
    tleg = arrayfun(@(t) sprintf('t=%s %s',formatsci(t/tunit,'eco'),tunitleg),tfit,'UniformOutput',false);
    
    % build axes
    if issample % basic plot
        wholedomain = sample;
        positionvar = 'x';
        positionlabel = 'distance';
        legposition = 1;
        hs = [NaN;subplots(1,[1 1],.04,.1)];
    else % as above + skin in position hs(1)
        wholedomain = 'apple';
        positionvar = 'r';
        positionlabel = 'radial position';
        legposition = 2;
        hs = [
            subplots([1 1],[1 1],.04,.1,'alive',[1 3])
            subplots(1,[1 1],.04,.1,'position',gcf,'alive',2)
            ];
        % concentration profiles through the skin
        subplot(hs(1))
        hp = plot(Rout.skin.x/lunit,interp1(Rout.skin.t,Rout.skin.C',tfit,'pchip')','-','linewidth',1.5);
        arrayfun(@(i) set(hp(i),'color',col(i,:)),1:o.plotoptions.nprofiles)
        hl = legend(hp(1:nleg),tleg(1:nleg),1,'fontsize',8); set(hl,'box','off')
        xlabel(sprintf('distance to surface (%s)',o.plotoptions.lunit),'fontsize',16)
        ylabel('concentration (?)','fontsize',16)
        title({'Distribution of the contaminant' 'through the skin'},'fontsize',10,'fontweight','bold')
    end
    
    % concentration profiles in the apple
    subplot(hs(2))
    hp = plot(Rout.(wholedomain).(positionvar)/lunit,interp1(Rout.(wholedomain).t,Rout.(wholedomain).C',tfit,'pchip')','-','linewidth',1.5);
    arrayfun(@(i) set(hp(i),'color',col(i,:)),1:o.plotoptions.nprofiles)
    hl = legend(hp(1:nleg),tleg(1:nleg),legposition,'fontsize',8); set(hl,'box','off')
    xlabel(sprintf('%s (%s)',positionlabel,o.plotoptions.lunit),'fontsize',16)
    title(sprintf('Profiles trough the whole %s',wholedomain),'fontsize',10,'fontweight','bold')
    
    % concentration kinetics
    subplot(hs(3)), hold on
    if issample % one kinetic with sample
        plot(Rout.(sample).t/tunit,Rout.(sample).Cav,'-','linewidth',1.5,'color',rgb(o.plotoptions.color.(sample)))
        ylabel(sprintf('concentration in the %s (?)',sample),'fontsize',16)
    else % three kinetics with the full apple
        hp = [
            plot(Rout.skin.t/tunit,Rout.skin.Cav,'-','linewidth',1.5,'color',rgb(o.plotoptions.color.skin))
            plot(Rout.skin.t/tunit,Rout.flesh.Cav,'-','linewidth',1.5,'color',rgb(o.plotoptions.color.flesh))
            plot(Rout.skin.t/tunit,Rout.apple.Cav,'-','linewidth',2,'color',rgb(o.plotoptions.color.apple))
            ];
        hl = legend(hp,{'skin' 'flesh' 'apple'},2,'fontsize',12); set(hl,'box','off')
        ylabel('concentration (?)','fontsize',16)
    end
    xlabel(sprintf('time (%s)',o.plotoptions.tunit),'fontsize',16)
    
    % format axes
    formatax(hs(~isnan(hs)),'fontsize',12)
    
%     if ~noverbosity
%         display comments for printing (note that filename is a property assigned when the figure was created with formatfig)
%         dispf('\nTo print the generated figure, use one the command (pwd=local directory)')
%         dispf( '\t>print in PDF\t   print_pdf(600,get(gcf,''filename''),pwd,''nocheck'')')
%         dispf( '\t>print in PNG\t   print_png(400,get(gcf,''filename''),pwd,'''',0,0,0)')
%         dispf('\t>export as CSV\t   print_csv(get(gcf,''filename''),pwd)')
%         dispf('\t>export as FIG\t   saveas(gcf,fullfile(pwd,sprintf(''%%s.fig'',get(gcf,''filename''))))')
%     end
    
end