function S=senspatankar_wrapper(S0,varargin)
%SENSPATANKAR_WRAPPER wrapper of senspatankar to use all inputs with real units instead of dimenesionless times and lengths
%  SYNTAX: Sdimensionless=senspatankar_wrapper(Swithunits [,'property1',value1,...,keyword,...])          
%          => to be used in R = senspatankar(Sdimensionless) or in any similar function (see section "see also")
%  INPUTS:: properties to be rescaled
%           S0.t or 't': time scale with units in s
%           S0.L or 'L': dilution ratio related to the whole thickness (VP/VF)
%  INPUTS:: Specific property/value and keyword
%           'iref' enables to force a particular layer as reference
%           'tunit' time unit for conversion (default = 'days')
%           'lunit' length unit for conversion (default = 'um')
%          'tlabel' time label (default = 'time')
%          'llabel' position label (default = 'position')
%  INPUTS:: keywords
%           'display' displays the reference layer and its related conductance conductance
%  INPUTS:: properties of senspatankar (to be overdefined)
%       example: senspatankar_wrapper(S0,'k',[10 1 1 1])
%  Note: parameters defined after S0 have higher precedence
%       (example: if S0.options exists and 'options' is defined, the value of 'options' is used instead of S0.options)
%
%  OUTPUTs: Sdimensionless has the same fields as S0 (+ user overrides)
%       Modified fields
%                t: dimensionless Fourrier time related to l(iref) and D(iref)
%                L: rescaled dilution factor related to l(iref)
%       Added fields
%             iref: index of the reference layer with the lowest conductance
%        timescale: scale to convert t in s
%      lengthscale: scale to convert x in m
%   hasbeenwrapped: true (flag)
%      wrap.time and wrap.position with fields (see EXAMPLE for details)
%           scale: anonymous function @(R) giving the scale matching unit
%                  where R = senspatankar(Sdimensionless) or the output of any similar function
%            unit: string giving unit
%           label: string giving label
%   
%
%   see also: SENSPATANKAR, SENSPATANKARC, SETOFFPATANKAR, SENSPATANKAR_RESTART, SENSPATANKARCGEOMETRY
%
%
%   EXAMPLE (note that FMECAunit is used to enter real units)
%{
    F	= 	struct(...
                    'Bi'		, 	1e3,...	Biot [hm.l(iref)/D]
                    'k'			,	[5 1],... ki, i=1 (layer in contact with the liquid)
                    'D'         ,   [3e-15 1e-14],... diffusion coefficient (m2/s)
                    'k0'        ,   1,... 0 = liquid
                    'l'         ,   FMECAunit('l',{'0.1 mm' '40 um'}),...
                    'L'			,	1/50,...	dilution factor (respectively to iref)
                    'C0'		,	[0 500],...	initial concentration in each layer
                    't'         ,   FMECAunit('t',{'1 h' '10 h' '1 d' '1 w' '1 mo' '1 y'}) ...
                        );
    S = senspatankar_wrapper(F,'display','lunit','mm','tunit','months')
    R = senspatankar(S);
    subplot(121), plot(R.t*S.wrap.time.scale(R),R.CF), xlabel(S.wrap.time.label), ylabel('Concentration (?)')
    subplot(122), plot(R.x*S.wrap.position.scale(R),R.Cx'), xlabel(S.wrap.position.label)
%}

% MS-MATLAB-WEB 1.0 - 13/05/09 - Olivier Vitrac - rev. 28/11/2015

% revision history
% 28/04/11 add S.lengthscale
% 23/10/15 improved help
% 24/10/2015 add hasbeenwrapped (i.e. dimensionless)
% 27/11/2015 major update, user override for iref and all properties of senspatankar
% 28/11/2015 extract options correctly

%% default
keywords = 'display';
default = struct('iref',[],'options',[],'tunit','days','lunit','um','tlabel','time','llabel','position');
usertscale = 24*3600; % it should match default.tunit
userlscale = 1e-6;    % it should match default.lunit
mandatoryfields = {'t' 'L' 'l' 'k' 'D'};
isvector = [false false true true true];

%% arg check
% 1) extract the specific options and keywords of the function
% note that options (for ODE) are extracted at this stage with 'nostructexpand'
[options,remain] = argcheck(varargin,default,keywords,'nostructexpand');
if length(options.iref)>1, error('iref is a scalar property'); end
missingfields = setdiff(mandatoryfields,fieldnames(S0));
if ~isempty(missingfields)
    dispf('\n======== ERROR SENSPATANKAR_WRAPPER ========',length(missingfields)),
    cellfun(@(f) dispf('\t missing property: %s',f),missingfields)
    error('SENSPATANKAR_WRAPPER:: %d properties are missing',length(missingfields)),
end
% 2) populate S with overdefined user parameters
S = argcheck(remain,S0,'','keep');
if ~isempty(options.options), S.options = options.options; end
n = unique(cellfun(@(f) length(S.(f)),mandatoryfields(isvector)));
if length(n)>1, error('the size of l, k and D are not consistent'); end

%% main
[crit,irefcandidate] = min(S.D./(S.k.*S.l));  % PERMEATION CRITERION
if ~isempty(options.iref) && ~isnan(options.iref)
    if options.iref>0 && options.iref<=n
        iref = options.iref;
    else
        error('iref should be an integer ranging between %d and %d',1,n)
    end
else
    iref = irefcandidate;
end
S.iref = iref;
S.L = S.L * S.l(iref) / sum(S.l);
S.t = S.t*S.D(iref)/S.l(iref)^2;
S.timescale = S.l(iref)^2/S.D(iref);
S.lengthscale = S.l(iref); % reference length scale (added 28/04/11)
S.hasbeenwrapped = true;

%% result wrapper
if ~strcmp(options.tunit,default.tunit) % update units if different units are requested
    usertscale = FMECAunit('t',sprintf('1 %s',options.tunit));
end
if ~strcmp(options.lunit,default.lunit) % update units if different units are requested
    userlscale = FMECAunit('l',sprintf('1 %s',options.lunit));
end
usertlabel = sprintf('%s (%s)',options.tlabel,options.tunit);
userllabel = regexprep(sprintf('%s (%s)',options.llabel,options.lunit),'um','\\mum');
S.wrap = struct(...
    'time',struct('scale',@(R)R.timebase/usertscale,'unit',options.tunit,'label',usertlabel),...
    'position',struct('scale',@(R)R.F.lengthscale/userlscale,'unit',options.tunit,'label',userllabel) ...
    );

%% diplay
if options.display
    dispf('\n%s::\nthe reference layer is the %dth of %d (conductance %0.4g m/s)',mfilename,iref,n,crit)
end