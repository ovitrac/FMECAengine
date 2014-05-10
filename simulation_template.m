% Example to launch a batch of simulations using SFPP3 engines
% INRA\Olivier Vitrac - rev. 08/05/14
%
% Revision history
% 26/04/14 Release candidate (based on SFPP3 engines only)
% 05/05/14 add thermodynamical data, add experimental data
% 08/05/14 It includes the new FMECAengine version v0.5


clc, close all
%% General definitions
% Private definitions (object-oriented programming)
engine  = @(simobject)    catstruct(senspatankar(senspatankar_wrapper(simobject)),struct('Fo',[],'nfo',''));
setprop = @(prop,default) argcheck(prop,default,'','case','keep');
default = senspatankar;

% User section (set path, outputfolder, figure name for printing)
switch localname % add your machine
    case 'WSLP-OLIVIER3', root = 'C:\Data\Olivier\INRA\Projects\SaintGobain2013'; local = fullfile(root,'Data','hafid');
    case 'lp-mol4.agroparistech.fr', root = '/home/olivier/hafid/SaintGobain2013'; local = fullfile(root,'matlab');
    case 'lp-mol3.agroparistech.fr', root = '/home/olivier/gwenael/SaintGobain2013'; local = fullfile(root,'matlab');
    case 'your machine', local = 'your full path';
    otherwise, error('add a line case ''your machine'', local = ''your full path'';')
end
thermodatafile = fullfile(root,'data','sample_updated.ods');
experimentaldatafile = fullfile(root,'data','hafid','march2014','resultats.ods');    
outputfolder = 'simulations'; if ~exist(fullfile(local,outputfolder),'dir'), mkdir(local,outputfolder), end
figname = sprintf('batchsimulations_%s',datestr(now,'dd-mmm-yy_HH:MM'));
paperposition = [0.1920    9.4226   20.6000   10.8323]; % paper position (use print preview, retrieve values with get(gcf,'paperposition))
position      = [179  371   1424  609];                 % window position
PRINTON = false; % set it to true to print results as PDF and PNG

% fmecaengine definitions (for HTML outputs)
enginefmeca = @(fmecaobject) fmecaengine('fmecamainfile',fmecaobject,'nograph',false,'noprint',false,'mergeoutput',true,'cls',false,...
                            'local',local,'outputpath',fullfile(outputfolder,'fmecaoutput'));

%% Thermodynamical properties
R = 8.3145; % J/mol/K
% v02(T in K, p in Pa) -- PVT relationships and equations of state of polymers (Eq. of Tait)
V0T = @(T) 0.7196 + 5.581e-5 * T + 1.468e-6 * T.^2; % cm3/g with T in degC
B0T = @(T) 3290.83 * exp(-2.321e-3-5.321e-3*T); % MPa with T in degC
C   = 0.0894; % (Cutler's constant)
p   = 0.1; % MPa
v02 = @(T) V0T(T) * (1 - C * log(1+p./(C*B0T(T)))); % cm3/g
% solutes
solprop = loadodsprefetch(thermodatafile,'sheetname','prop','prefetchprefix','PREFECTCH_properties_');
solpropA = bykeywords(solprop,'plasticizer','plasticizer A');
solpropB = bykeywords(solprop,'plasticizer','plasticizer B');
v01A = @(T) interp1(solpropA.tempK,solpropA.molarvolume,T+273)/solpropA.Mw(1); % cm3/g with T in degC
v01B = @(T) interp1(solpropB.tempK,solpropB.molarvolume,T+273)/solpropB.Mw(1); % cm3/g
% volume fraction in polymer
% w1Aw2 = mass A / mass dry P
% w1Bw2 = mass B / mass dry P
phi2 = @(T,w1Aw2,w1Bw2) v02(T) ./ ( v02(T) + w1Aw2.*v01A(T) + w1Bw2.*v01B(T) );
% apparent density of the polymer
dapp2 = @(T,w1Aw2,w1Bw2) phi2(T,w1Aw2,w1Bw2) ./ v02(T);
% apparent density of ethanol
dEtOH = @(T) 0.8*exp(-(T-25)*750e-6);
% conversion of K = CF/CP
% from: mass reference frame = dry mass of P
% to: volume concentrations
Kdrymass2Kvol = @(K,T,w1Aw2,w1Bw2) K .*  dEtOH(T) ./ dapp2(T,w1Aw2,w1Bw2);


%% Reference simulation (note that not explicit parameters are inherited from default)
% composition data
plasticizercontent = 0.45;
CP0      = plasticizercontent/(1-plasticizercontent); % total amount of plasticizer in polymer
CP0ratio = [0.85 0.15]; 
% main properties are explicitely defined (others are inherited from default)
tscale = @(t) logspace(-3,log10(t),1e4)'; hours = 3600;
refsim = default;
refsim.Bi = 1e9;
refsim.C0 = CP0ratio(1)*CP0;     % initial concentration
refsim.D  = 1.2e-12; % m2/s
refsim.k  = Kdrymass2Kvol(0.3,25,CP0ratio(1)*CP0,CP0ratio(2)*CP0);   % KF/P (when k0=1)
refsim.L  = 1/6;   % LP/F (volumic)
refsim.l  = 6.35e-3; % m
refsim.t  = tscale(96*hours);% s (logscale)
refsim.nfo = 'default simulation Xiaoyi';

% Private definitions to be used with FMECAengine (it requires FMECAengine version 0.50 or above)
% Note that we use shorthands instead of the conventional fields of FMECAengine
fmecadefault = @(refsim,relthick,t1,id) fmecaengine('fmecamainfile',...
{   'idusercode',{sprintf('%s01',id) sprintf('%s02',id)},...
    'nlayers',2,'nsteps',2,'constructor','l',repmat({refsim.l*[relthick (1-relthick)]},1,2),... Two steps x Two layers
    'D',{refsim.D*[1 1] refsim.D*[1/100 1]},... Two steps x Two layers
    'KFP',repmat({refsim.k*[1 1]},1,2),... Two steps x Two layers
    'CP',repmat({refsim.C0*[1 1]},1,2),... Two steps x Two layers
    'ts',[t1 refsim.t(end)],... Two steps
    'parent',{'' sprintf('%s01',id)},... Two steps: chaining rumles (steps are named STEP01 STEP02)
    'lFm',refsim.l/refsim.L,'Binounit',refsim.Bi,...
    'relthick',relthick,...
    't1',t1 ...
    });

%% build batch and launch simulations when needed
% TIP: use  setprop({'property1',value1, 'property2',value2,...},refsim) to change properties of refsim
% ========================= S I M U L A T I O N   D E F I N I T I O N S ============================================
batch = [
    refsim
    %setprop({'k',Kdrymass2Kvol(0.036,25,CP0ratio(1)*CP0,CP0ratio(2)*CP0),'D',2.88e-10,'t',tscale(90*hours),'nfo','hafid'},refsim)
    setprop({'k',Kdrymass2Kvol(0.0136/3,25,CP0ratio(1)*CP0,CP0ratio(2)*CP0),'D',8.2e-10,'t',tscale(90*hours),'nfo','KA=0.0136 configuration cylindrique'},refsim)
    %setprop({'k',Kdrymass2Kvol(0.0136,25,CP0ratio(1)*CP0,CP0ratio(2)*CP0),'D',8.2e-10,'t',tscale(90*hours),'L',refsim.L*2,'nfo','L=1/3'},refsim)
    %setprop({'k',Kdrymass2Kvol(0.0136/2,25,CP0ratio(1)*CP0,CP0ratio(2)*CP0),'D',8.2e-10,'t',tscale(90*hours),'nfo','k/2'},refsim)
    %setprop({'k',Kdrymass2Kvol(0.075,25,CP0ratio(1)*CP0,CP0ratio(2)*CP0),'D',5.63e-9,'t',tscale(90*hours),'nfo','gwenael'},refsim)
    %setprop({'k', 0.05,'nfo',',k=0.05'},refsim)
    %setprop({'k', 0.01,'nfo','k=0.01'},refsim) % modify 1 parameter
    %setprop({'l',[1e-7 1]*refsim.l,'D', [1e-16 1e-12],'k',[.3 .03],'C0',[1 1],'nfo',sprintf('barrier layer\nrelative thick.10^{-7}')},refsim)
    %setprop({'l',[1e-6 1]*refsim.l,'D', [1e-16 1e-12],'k',[.3 .03],'C0',[1 1],'nfo',sprintf('barrier layer\nrelative thick.10^{-6}')},refsim)
    %setprop({'l',[1e-5 1]*refsim.l,'D', [1e-16 1e-12],'k',[.3 .03],'C0',[1 1],'nfo',sprintf('barrier layer\nrelative thick.10^{-5}')},refsim)
    ]; % add simulations here
% ==================================================================================================================

% ========================== launch simulations (only for thoses that need to be updated) ==========================
nbatch = length(batch);
t0 = clock;
if ~exist('savebatch','var') || ~exist('res','var') || isempty(res); savebatch = []; res=struct([]); end
dispf('\nBATCH of %d SIMULATIONS (%s)\n%s',nbatch,datestr(now),repmat('-',1,45))
nupdated = 0;
for ibatch = 1:nbatch
    t1 = clock;
    if isempty(savebatch) || ibatch>length(savebatch) || ~structcmp(batch(ibatch),savebatch(ibatch)) || ...
            isempty(res) || ~isstruct(res) || ~isfield(res,'C') || length(res)<nbatch || isempty(res(ibatch).C)
        if isempty(batch(ibatch)), res(ibatch).nfo = sprintf('batch %d',ibatch); else res(ibatch).nfo = batch(ibatch).nfo; end
        dispf('BATCH %d/%d running...%40s',ibatch,nbatch,regexprep(batch(ibatch).nfo,'\n',' '))
        if ibatch==1
            res = engine(batch(1)); 
            res.Fo = res.t; % Fourier times
            res.t  = res.Fo * res.timebase; % real times
        else
            res(ibatch) = engine(batch(ibatch)); %#ok<*SAGROW>
            res(ibatch).Fo = res(ibatch).t; % Fourier times
            res(ibatch).t  = res(ibatch).Fo * res(ibatch).timebase; % real times
        end
        if isempty(batch(ibatch)), res(ibatch).nfo = sprintf('batch %d',ibatch); else res(ibatch).nfo = batch(ibatch).nfo; end
        nupdated = nupdated+1; dt = etime(clock,t0);
        dispf('\tend in %0.4g s (elapsed time %0.4g s, remaining time %0.4g s)',etime(clock,t1),dt,dt*(nbatch/ibatch-1))
    else
        dispf('nothing to do for BATCH %d/%d (%s)',ibatch,nbatch,res(ibatch).nfo)
    end
end
dispf('%d/%d simulations updated in %0.5g s',nupdated,nbatch,etime(clock,t0))
savebatch = batch;

%% ========================== add experimental data ==========================
% load data and data collection
% Data model: data(itemperature).raw(ikinetic).field with field = t, CFA, CFB
raw = loadodsprefetch(experimentaldatafile,'sheetname','all');
data = struct(...
    'temperature',{25 40},... temperature
    'title',{'Contact 25^\circC' 'Contact 40^\circC'},... corresponding titles
    'kinetics',{{'Gwenael25' 'Hafid25' 'kineticswithtwodifferentvolume'} {'Gwenael40'}},... sheet names
    'marker',{{'x' '+' '*'}},...
    'markerapprox',{{'o' 's' '^'}},...
    'legend',{{'bottom' 'wall' 'wall'} {'bottom'}},... corresponding legends
    'hmin',{[1000 1000 1000] [4000]},...
    'hwidth',{[5 5 5] [20]} ...
    );
% Populate raw fields
ndata = length(data);
for i=1:length(data)
    data(i).nkinetics = length(data(i).kinetics);
    for j=1:data(i).nkinetics
        % raw
        data(i).raw(j) = raw.(data(i).kinetics{j});
        rawfields = fieldnames(data(i).raw)';
        valid = find(~isnan(data(i).raw(j).mAm0A));
        [~,isort] = sort(data(i).raw(j).t(valid));
        nozero = ~any(data(i).raw(j).t==0);
        % add zeros if needed, sort data
        for f=rawfields
            data(i).raw(j).(f{1}) = data(i).raw(j).(f{1})(valid(isort));
            if nozero
                if isnumeric(data(i).raw(j).(f{1})), data(i).raw(j).(f{1}) = [0;0;0;data(i).raw(j).(f{1})];
                else data(i).raw(j).(f{1}) = [{'zero';'zero';'zero'};data(i).raw(j).(f{1})];
                end
            end
        end % next f
    end % next j
end % next i

%% ========================== generalized simulations ==========================
t1 = 7*hours;
fmecasim = [fmecadefault(refsim,0.001,t1,'low');fmecadefault(refsim,0.01,t1,'high')];
fmecares = enginefmeca(fmecasim);


%% ========================== do plots ==========================
currentfig = gcf; clf; formatfig(gcf,'figname',figname,'paperposition',paperposition,'position',position)
hs = subplots([1 1.3],1,0); hp = zeros(nbatch,1); leg = cell(nbatch,1); colors = cbrewer('qual','Set2',nbatch);
for i=1:nbatch
    L = res(i).F.L*res(i).F.lrefc(end);
    CPO = sum(res(i).F.C0.*res(i).F.lref)/res(i).F.lrefc(end);
    CF_LCP0 = res(i).CF/(L*CP0);
    subplot(hs(1)), hold on, plot(res(i).t/hours,CF_LCP0,'-','color',colors(i,:),'linewidth',2);
    subplot(hs(2)), hold on, hp(i) = plot(sqrt(res(i).t/hours),CF_LCP0,'-','color',colors(i,:),'linewidth',2);
end
% add generalized simulations
fmecacolors = rgb('bluecolors');
leg = {res.nfo};
for i=1:length(fmecares)
    CF_LCP0 = fmecares(i).CF/(L*CP0);
    subplot(hs(1)), hold on, plot(fmecares(i).t/hours,CF_LCP0,'-','color',rgb(fmecacolors{i}),'linewidth',2);
    subplot(hs(2)), hold on, hp(end+1) = plot(sqrt(fmecares(i).t/hours),CF_LCP0,'-','color',rgb(fmecacolors{i}),'linewidth',2);
    leg{end+1} = regexprep(sprintf('+%s \n',fmecares(i).path{:}),'\n$','');
end
% legend and labels
hl = legend(hp,leg,'location','NorthEastOutside','fontsize',10); set(hl,'box','off')
xticklabels = cellstr(get(hs(1),'xticklabel')); xticklabels{end} = ''; set(hs(1),'xticklabel',xticklabels)
subplot(hs(1)), xlabel('time (hours)','fontsize',14), ylabel('C_F/(L\cdotC_P)','fontsize',16)
xax = xlim;
plot(data(1).raw(2).t/hours,data(1).raw(2).CFA*dEtOH(25)/L*CPO,'r+') %hafid results at 25°C for A
%plot(data(1).raw(2).t,data(1).raw(2).mBm0B,'r+') %hafid results at 25°C for B

set(hs(1),'xlim',xax)
subplot(hs(2)), xlabel('sqrt(time) (s^{1/2})','fontsize',14)
formatax(hs,'fontsize',10), titles(hs,'','suffix','.','x',.1,'y',.9,'fontsize',8,'fontweight','bold')


%% ========================== print ==========================
if nupdated>0 && PRINTON
    print_pdf(600,get(gcf,'filename'),fullfile(local,outputfolder),'nocheck')
    print_png(300,get(gcf,'filename'),fullfile(local,outputfolder),'',0,0,0)
end

