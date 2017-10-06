function [R, profile, Sout] = tube4fmecaengine(S0,varargin)
% tube4fmecaengine interprets additional arguments for opened system, transforms into fmecaengine inputs and launch simulation
% Syntax:
% R = tube4fmecaengine(S0,'property1',value1,'property2',value2);
%
%               S0 : scalar structure or 1xnsections structure array coding for the materials properties of each section
%                    default constructor:
%                           S0 = fmecaengine2('fmecamainfile',{'constructor'})
%                           S0 = fmecaengine2('fmecamainfile',{'constructor','parameter1',value1,'parameter2',value2,...})
%                    For simple simulations, S0 can be defined as 
%                           S0.parameter1 = value1; S0.parameter2 = value2; etc.
%                           S0={'parameter1',value1,'parameter2',value2,...}
%                   To be implemented: SO is a valid ODS file as FMECAengine2() can manage (+ add the right sheet)
%                   NOTE1: contact times (ts), liquid thickness (lFm) will be calculated by tube4fmecaengine according to flowrate, sectionlength and crosssectionarea
%                   NOTE2: S0 can code for sections with different number of layers (use NaN for non-used properties)
%                   NOTE3: if length(S0) does not match the number of sections (see sectionlength), missing values are padded by repeating the last definition (argpad)
%
%   Pair property/value coding for tube properties, flow and related physical approximations
%          flowrate: flow rate in SI unit (scalar value, default = 1e-8 m3/s), STEADY STATE IS UNIT
%     sectionlength: 1xnsections array coding for the length of each section (default = 1 m)
%                    NOTE3: length(sectionlength) determines the number of sections
%  crosssectionarea: scalar or 1xnsections array (default = 1e-6 m2), missing values are padded by repeating the last value (argpad)
%             ntmax: scalar value, numbers division in time scale, ntmax must be, at least, equal to length(sectionlength) (default = [])
%
%   Pair property/value coding for project
%           session: session name of pcurrent project (default = autoprojectname(5,true)) 
%            ploton: plot kinetics if true (default = true)
%
%  OUTPUT
% R: output structure of fmecaengine2
% PROFILE: structure with fields
%          t: 1 x ntmax array of contact time (secondes)
%         CF: as CF output of R, 1 x ncirculation array of conc. in fluid in each section at all time (kg/m3)
%         mF: 1 x ntmax array of mass in each section at all time (kg)
%        CPi: 1 x ntmax array of conc. in polymer at t end  (kg/m3)
%        CFi: 1 x ntmax array of conc. in fluid at t end in each section (kg/m3)
% ex. R = tube4fmecaengine
%{
      S0(1) = struct('nlayers',2,'nsteps',1,'idusercode',{'step1'},'l',[6e-4 NaN],'D',[1e-11 NaN],'KFP', 1,'CP',[525 NaN],'CF0',0,'ts',1000,'lFm',6.7e-4,'Binounit',50);
      S0(2) = struct('nlayers',2,'nsteps',1,'idusercode',{'step1'},'l',[6e-4 1e-3] ,'D',[1e-14 1e-12],'KFP', 0.1,'CP',[0 400],'CF0',0,'ts',1000,'lFm',6.7e-4,'Binounit',50);
      project = 'armed'; 
      [R1, profile1, S1] = tube4fmecaengine(S0,'flowrate',20e-6/3600,'sectionlength',[0.1193    0.2330    0.3561    0.4792    0.6024],...
                                            'session',project,'crosssectionarea', 6e-6,'ntmax',100,'ploton',true);
%}
%
% UPDATE May 2017 (to be used for debugging)
%{
S = struct('nlayers',1,'nsteps',1,'idusercode','step1','l',1e-3,'D',1e-14,'KFP', 0.1,'CP',1000,'CF0',0,'ts',50,'lFm',6e-3,'Binounit',50);
F1 = struct('session','flowrate1','flowrate',1e-8,'length',0.1,'nsections',2,'crosssectionarea',pi*1e-3^2);
[R1, profile1, S1] = tube4fmecaengine(S,F1,'ntmax',100,'ploton',true);
%}


% Migration 2.1 - 22/02/2016 - INRA\Olivier Vitrac - Mai Nguyen - rev. 30/08/2017
% Revision history
% 07/03/2016 add 'savememory' 
% 08/03/2016 add CFi and CPi in output
% 10/05/2017 fix mF and add help
% 30/08/2017 add source step

%% DEBUGMODE
% use dbstop in fmecaengine2 at 1069
DEBUGMODE = false;

%% default
% Static properties
Sdefault = struct('nlayers',1,...
                 'nsteps',1,...
                 'idusercode',{'step1'},...
                 'l',1e-4,...
                 'D',1e-12,...
                 'KFP', 1,...
                 'CP',300,...
                 'CF0',0,...
                 'ts',300,...
                 'lFm',1e-4,...
                 'Binounit',100);
% Dynamic properties
Fdefault = struct('flowrate',1e-8,...
                  'length',.1,...
                  'nsections',3,...
                  'crosssectionarea',1e-6,...
                  'ntmax',100,...
                  'session',autoprojectname(5,true),...
                  'ploton',true); % in SI units

%% argcheck
% dynamic inputs
oF = argcheck(varargin,Fdefault);
% static inputs
oS = argcheck(S0,Sdefault);
if nargin<1, S0 = []; end
if isempty(S0), S0 = Sdefault; end
if length(S0) ~= oF.nsections, S0 = argpad(S0,oF.nsections); end
% additional control
sectionlength = oF.length/oF.nsections;
vsection = sectionlength .* oF.crosssectionarea;
dt = vsection/oF.flowrate; % residence time in each section
ntmin = oF.nsections;
nt = max(ntmin,min(oF.ntmax,ceil(oS.ts/dt)));
% building of fmecaengine constructor
template = fmecaengine2('fmecamainfile',{'constructor','session',oF.session,S0(1),'parent','','parentF',''});
S = repmat(template,oF.nsections*nt - oF.nsections*(oF.nsections-1)/2,1);
source = argcheck(struct('idusercode','source','CF0',oS.CF0,'isstatic',true),template,'','case');
stepname = @(it,iz) sprintf('t%04dz%04d',it,iz);
istep=0;
for it = 1:nt
    istep = istep + 1;
    % section in contact with the source (inlet)
    if it>1, parentinlet = stepname(it-1,1); else, parentinlet = ''; end
    S(istep) = fmecaengine2('fmecamainfile',{'constructor','session',oF.session,S0(1),'parent',parentinlet,'parentF','source'});
    S(istep).idusercode = stepname(it,1);
    S(istep).ts = dt;
    % sections beyond the first
    for iz = 2:oF.nsections
        if it>=iz
            istep = istep + 1;
            S(istep) = fmecaengine2('fmecamainfile',{'constructor','session',oF.session,S0(iz),'parent','','parentF',''});
            S(istep).idusercode = stepname(it,iz);
            S(istep).ts = dt;
            if it>iz 
                S(istep).parent = stepname(it-1,iz);
            else
                S(istep).parent = ''; % no parent before the first step
            end
            if iz>1
                S(istep).parentF = stepname(it-1,iz-1);
            else
                S(istep).parentF = ''; % no parentF for the first
            end
        end
    end
end
S = [source;S];

% launch simulations
t0 = clock;
if DEBUGMODE
    R = fmecaengine2('fmecamainfile',S(:),'ramdisk','nograph',true,'noprint',true,'nohtml',true); %'noplot','savememory');
else
    R = fmecaengine2('fmecamainfile',S(:),'ramdisk','nograph',true,'noprint',true,'nohtml',true,'noplot','savememory');
end
dispf('\n>> Total execution time: %0.4g s',etime(clock,t0))

% create profile in each section
profile = repmat(struct('nt',nt,'t',[],'CF',[],'CFcum',[],'CFi',[],'CPi',[],'mF',[]),oF.nsections,1);
for iz = 1:oF.nsections
    profile(iz).t  = (iz:nt)*R.(stepname(iz,iz)).t; % s residence time * ntmax    R.(stepname(iz,iz)).t = dtsub(iz)
    profile(iz).CF = arrayfun(@(it) R.(stepname(it,iz)).CF,iz:nt); % kg/m3 
    profile(iz).CFi = arrayfun(@(it) R.(stepname(it,iz)).CFi,iz:nt); % kg/m3
    profile(iz).CPi = arrayfun(@(it) R.(stepname(it,iz)).CPi,iz:nt); % kg/m3
    profile(iz).mF = cumsum(profile(iz).CF*vsection); % kg cumtrapz(profile(iz).t,profile(iz).CF*o.flowrate);
    profile(iz).CFcum = profile(iz).mF ./ (profile(iz).t*oF.flowrate);
end

% plot
if oF.ploton
    figure, hs = subplots(1,[1 1],0,0);
    subplot(hs(1)), plotpub({profile.t},{profile.mF}), ylabel('m_F (kg)')
    subplot(hs(2)), plotpub({profile.t},{profile.CF}), ylabel('C_F (kg\cdotm^{-3})') , xlabel('time (s)')
    formatax(hs)
end
if nargout > 2, Sout = S; end
