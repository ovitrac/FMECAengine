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
%           'display' displays the reference layer and its related  conductance
%  INPUTS:: pair property/value for execution and evaluation (default=[])
%           'run' sets an anonymous function to run output S
%           'eval' sets an anonymous function to evaluate the previous result at some particular conditions (time and position)
%           TIP: uses 'run' and 'eval' for fitting (see advanced example)
%  INPUTS:: properties of senspatankar (to be overdefined)
%       example: senspatankar_wrapper(S0,'k',[10 1 1 1])
%  Note: parameters defined after S0 have higher precedence
%       (example: if S0.options exists and 'options' is defined, the value of 'options' is used instead of S0.options)
%  Note: iref is an index based on existing layers (i.e. if C0 = [-1 100 0] with layer 1 removed, iref must be ranged between 1 and 2)
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
%
%
% ADVANCED EXAMPLE including fitting
%{
    expkin = struct(... migration data of 4-nitrophenol from the master of Sofiane (2016)
        't', [0    0.0153    0.0358    0.9589    3.9951    4.9951    5.9469    7.0947    8.0947   13.0551   26.0214   46.9040]',... days
       'CF', [0     0.1276    0.0947    0.3050    1.2825    1.6613    1.7807    2.0758    2.1898    2.8587    2.4089    2.8855] );
    inthefit = true(size(expkin.t));    
    days = FMECAunit('t','1 day');
    S0 = struct('Bi',1e6,...
                'k',1,...
                'D',1e-14,...
                'k0',1,...
                'l',150e-06,...
                'L',0.0783,...
                'C0',37.5,...
                'CF0',0,...
                't',(linspace(0,sqrt(50),1e4).^2)*days,...
                'run',@senspatankar ...
            );
    A0 = log10([S0.D S0.k]);
    evalR = @(res,inthefit) interp1(res.timebase*res.t,res.CF,expkin.t(inthefit)*days,'pchip');
    Rsim = @(A,evalR) senspatankar_wrapper(S0,'D',10^A(1),'k',10^A(2),'eval',evalR);
    crit = @(A,inthefit) norm(expkin.CF(inthefit)-Rsim(A,@(R) evalR(R,inthefit))).^2;
    fitoptions = optimset('display','iter','FunValCheck','on','MaxIter',1e3,'TolFun',1e-4,'TolX',1e-4);
    [Aopt,err] = fminsearch(@(A) crit(A,inthefit),A0,fitoptions);
    Ropt = Rsim(Aopt,[]);
    hp = plot(expkin.t,expkin.CF,'ro',Ropt.t*Ropt.timebase/days,Ropt.CF,'b-',expkin.t(inthefit),evalR(Ropt,inthefit),'bx');
    hl = legend(hp,{'experiment' 'fitted' 'calculated'},'location','southeast','fontsize',12); set(hl,'box','off')
    title(sprintf('\\rmD = \\bf%s\\rm m^2\\cdots^{-1}  -  K_{F/P} = \\bf%s',formatsci(10^Aopt(1)),formatsci(10^Aopt(2),'eco')),'fontsize',18)
    % as before, without end-1
    inthefit(end-1) = false;
    [Aopt2,err2] = fminsearch(@(A) crit(A,inthefit),A0,fitoptions);
    Ropt2 = Rsim(Aopt2,[]);
    hp2 = plot(Ropt2.t*Ropt2.timebase/days,Ropt2.CF,'g-',expkin.t(inthefit),evalR(Ropt2,inthefit),'g+');
%}


% MS-MATLAB-WEB 1.0 - 13/05/09 - Olivier Vitrac - rev. 14/04/2017

% revision history
% 28/04/11 add S.lengthscale
% 23/10/15 improved help
% 24/10/2015 add hasbeenwrapped (i.e. dimensionless)
% 27/11/2015 major update, user override for iref and all properties of senspatankar
% 28/11/2015 extract options correctly
% 30/11/2015 keep case on S0, iref based only on existing layers only
% 01/12/2015 fix the size control of C0
% 14/04/2017 add 'run' and 'eval'

%% default
keywords = 'display';
default = struct('iref',[],'options',[],'tunit','days','lunit','um','tlabel','time','llabel','position','run',[],'eval',[]);
usertscale = 24*3600; % it should match default.tunit
userlscale = 1e-6;    % it should match default.lunit
mandatoryfields = {'t' 'L' 'l' 'k' 'D' 'C0'};
isvector = [false false true true true true];

%% arg check
% 1) extract the specific options and keywords of the function
% note that options (for ODE) are extracted at this stage with 'nostructexpand'
[options,remain] = argcheck(varargin,default,keywords,'nostructexpand','case');
if length(options.iref)>1, error('iref is a scalar property'); end
missingfields = setdiff(mandatoryfields,fieldnames(S0));
if ~isempty(missingfields)
    dispf('\n======== ERROR SENSPATANKAR_WRAPPER ========',length(missingfields)),
    cellfun(@(f) dispf('\t missing property: %s',f),missingfields)
    error('SENSPATANKAR_WRAPPER:: %d properties are missing',length(missingfields)),
end
% 2) populate S with overdefined user parameters
S = argcheck(remain,S0,'','keep','case');
if isfield(S,'run') && isempty(options.run), options.run = S.run; end
if isfield(S,'eval') && isempty(options.eval), options.eval = S.eval; end
if ~isempty(options.options), S.options = options.options; end
n = unique(cellfun(@(f) length(S.(f)),mandatoryfields(isvector)));
if length(n)>1, error('the size of l, k D, C0 are not consistent (ranged between %d and %d)',min(n),max(n)); end

%% main
ok = find(S.C0>=0); nok = length(ok);
[crit,irefcandidate] = min(S.D(ok)./(S.k(ok).*S.l(ok)));  % PERMEATION CRITERION
if ~isempty(options.iref) && ~isnan(options.iref)
    if options.iref>0 && options.iref<=nok
        irefcandidate = options.iref;
    else
        error('iref should be an integer ranging between %d and %d (number of layers=%d)',1,nok,n)
    end
end
S.iref = irefcandidate;
iabsref = ok(irefcandidate);
S.L = S.L * S.l(iabsref) / sum(S.l); % normalization for all layers (even those removed with C0<0)
S.t = S.t*S.D(iabsref)/S.l(iabsref)^2;
S.timescale = S.l(iabsref)^2/S.D(iabsref);
S.lengthscale = S.l(iabsref); % reference length scale (added 28/04/11)
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
    dispf('\n%s::\nthe reference layer is the %dth of %d (conductance %0.4g m/s)',mfilename,iabsref,n,crit)
end

%% run if specified
if ~isempty(options.run)
     if isa(options.run,'function_handle')
        R = options.run(S);
    else
        error('eval must define an anonymous function')
     end
     % eval if required
    if ~isempty(options.eval)
        if isa(options.eval,'function_handle')
            S = options.eval(R);
            return
        else
            error('eval must define an anonymous function')
        end
    else
        S = R;
    end
end
