function res = senspatankar(F,ploton,dispon)
%SENSPATANKAR simulates transfer through n layers using a modified Patankar Method (see p 45)
%   the dimensionless formulation is similar to SENSN
%   all data are normalized according the reference layer or equivalently according to the layer with the lowest D/a value
%
%   IT IS THE RESPONSABILITY OF THE USER TO PROVIDE THE APPROPRIATE DIMENSIONLESS NUMBERS
%   >> a wrapper is used for the online version is available in ../www/home/diffusion_1DFVn.m
%   >> at command prompt, the following wrapper can be used
%           S = senspatankar_wrapper(Sreal) % see code of senspatankar_wrapper for details
%           R = senspatankar(senspatankar_wrapper(Sreal))
%
%   The full method is detailed in the following presentation: http://modmol.agroparistech.fr/sfpp3/SFPP3_Migresives/
%   Slides: 16-20
%       Finite volume method described in: http://en.wikipedia.org/wiki/Finite_volume_method      
%       Rankine-Hugoniot jump condition described in: http://en.wikipedia.org/wiki/Rankine%E2%80%93Hugoniot_conditions
%       Discretization scheme: http://modmol.agroparistech.fr/home/FVmultilayers.pdf 
%       Physics fully described in the following chapter: http://modmol.agroparistech.fr/pub/Vitrac_Hayert_Nova.pdf
%       Design of safe food packaging under uncertainty (in English) Authors: O. Vitrac and M. Hayert.
%       In "New trends chemical engineering research", Ed. L. P. Berton, Nova Science Publishers, New York, 2007, p251-292
%
%   Geometric convention
%       0=food or F, 1=layer in direct contact with food, 2, ..., n
%
%   Constructor: S = senspatankar; % returns the default input object
%     >with fields for a n layers structure (mandatory)
%          Bi: scalar, mass Biot number (relative to first layer: h*l1/D1)
%           k: 1xn vector, Henry-like coefficient for each layer (consequence: partition coefficient: K1/2 = C1eq/C2eq = k2/k1)
%           D: 1xn vector, diffusion coefficient for each layer (units in [l]^2/s)
%          k0: Henry-like coefficient for food (if k0=1, KF/i = ki with i=1..n)
%           l: 1xn vector, thickness of each layer (free unit)
%           L: LP/F * lref / sum(l), where LP/F = VP/VF (dimensionless number), VF=volume of food, VP=volume of the packaging
%          C0: 1xn, concentration in mass/volume for each layer
%           t: dimensionless times Fo=D*treal/lref^2 (Fourier number relative to lref) at which the solutions need to be calculated
%     >optional fields:default values
%    autotime:'on', optimize time steps for the time solver
%           N:1e4 , number of returned time steps (interpolated)
%       nmesh:200, number of finite element volumes
%    nmeshmin:10, minimum number of element volumes per layer
%     options:odeset('RelTol',1e-4,'AbsTol',1e-4,'Initialstep',1e-8,'Maxstep',.01,'Maxorder',2), structure defined by odeset
%        iref:force a different reference layer (automatically calculated field)
%      ivalid:index array of valid layers (default = true(1,n))
%
%   How to handle effects of density r0, r = [r1,r2,..ri]
%       keq = (k/k0) ./ (r/r0);
%   How to calculate the reference layer for S.L and S.t
%       equivalent conductance: a = D./keq
%       [~,iref] = min(a./k); lref = l(iref); Dref = D(iref);
%
%
%   To launch the simulation: R = senspantankar(S)
%   R = structure with fields
%          C: Nx1 vector, CP=f(t) kinetic in P (space averaged) at interpolated times R.t
%          t: Nx1 vector, times (dimensionless, use R.t*R.timescale for real units)
%          p: NaN (not used in senspatankar but for retrocompatibility with sensn)
%        peq: partial pressure at equilibrium
%             if k0=1, C0eq=peq
%             0eq = sum(S.lref*S.L.*S.C0)/(1+sum((S.k0./S.k .* S.lref*S.L)));
%          V: equivalent thickness of the packaging (reference scale = lref)
%         C0: S.nmeshx1 vector, initial solution
%          F: initial input S
%          x: 3*S.nmeshx1 vector, positions (control nodes and interface nodes) where the solution is calculated
%         Cx: Tx3*S.nmesh matrix, coding for C(t,x) where t is time and x is position, T=times at which the solution was really integrated
%         tC: Tx1 vector at which the solution was integrated (2nd order backward differences)
%         CF: Nx1 vector, CF=f(t) kinetic in food at interpolated times R.t
%         fc: Nx1 vector, cumulated fluw
%          f: Nx1 vector, flux
%     timebase: lref^2/Dref, use R.t*R.timebase to retrieve real time units
%
%   See also: senspatankar_wrapper, setoffpatankar, senspatankarT, senspatankarC, senspatankarnonlin (beta)

% MS-MATLAB-WEB 1.0 - 28/09/07 - Olivier Vitrac - rev 15/02/11

% Revision history
% 01/10/07 improve speed
% 16/03/09 add restart
% 17/11/10 improve documentation
% 15/02/11 fix Out of Memory for large Fo values (Fo>10)

% definitions
global timeout
timeout		= 800; % s
% options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Stats','yes','Initialstep',1e-5,'Maxstep',.05,'Maxorder',5);
options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Initialstep',1e-8,'Maxstep',.01,'Maxorder',2);
Fdefault	= 	struct(...
				'Bi'		, 	1e3,...	Biot [hm.L1/D]
				'k'			,	[1 1 1 1],...[0.5 3 2],...	ki, i=1 (layer in contact with the liquid)
                'D'         ,   [1e-16 1e-14 1e-14 1e-14],... diffusion coefficient
                'k0'        ,   1,... 0 = liquid
                'l'         ,   [50 20 10 120]*1e-6,...[50 20 10 120]*1e-6,... m
				'L'			,	200/1800,...	dilution factor (respectively to iref)
				'C0'		,	[0 500 500 500],...	initial concentration in each layer
				'options'	,	options...
					); % if iref is missing, it is indentified
Fdefault.t	= [0:.00001:.005 .01:.01:.1 .11:.1:5]; %0:.00001:.005; %[0:.00001:.005 .01:.01:.1 .11:.1:5]';
if Fdefault.Bi<10, Fdefault.t = Fdefault.t/(Fdefault.Bi/10); end
method		= 'cubic'; %'cubic';
ploton_default = false;
dispon_default = false;
nmesh_default  = 200; % number of nodes for a layer of normalized thickness 1
nmeshmin_default    = 10;
n_default	= 1e4;
MaxStepmax	= 10;
nstepchoice = 200;
zero = 1000*eps;

% arg check
initon = false;
if ~nargin, initon = true; end
if nargin<1, F = []; end
if nargin<2, ploton = []; end
if nargin<3, dispon = []; end
if isempty(F), F = Fdefault; end
if ~isfield(F,'autotime'), F.autotime = 'on'; end
if ~isfield(F,'n'), F.n = n_default; end
if ~isfield(F,'nmesh'), F.nmesh = nmesh_default; end
if ~isfield(F,'nmeshmin'), F.nmeshmin = nmeshmin_default; end
if ~isfield(F,'options'), F.options = options; end
if ~nargout, ploton=true; dispon=true; end
if isempty(ploton), ploton = ploton_default; end %#ok<NASGU>
if isempty(dispon), dispon = dispon_default; end

% physical check
m = Inf;
for prop = {'D' 'k' 'C0'};
    if strcmp(prop{1},'C0')
        F.(prop{1}) = F.(prop{1})(F.(prop{1})>=0);
    else
        F.(prop{1}) = F.(prop{1})(F.(prop{1})>0);
    end
    m = min(m,length(F.(prop{1})));
end
F.m = m;


if initon
    res = Fdefault;
    if dispon, disp('... init mode'), end
    return
end

if strcmpi(F.autotime,'on')
	ti		= linspace(min(F.t),max(F.t),F.n)';
else
	ti		= F.t;
end
MaxStepmax = max(MaxStepmax,ceil(F.t(end)/nstepchoice)); %#ok<NASGU>

% renormalization (fields X1->Xn)
F.k = F.k/F.k0; F.k0 = 1;   % by convention
a = F.D(1:m)./F.k(1:m);
if isfield(F,'iref')
    iref = F.iref;
else
    [crit,iref] = min( a./F.l(1:m) ); %#ok<ASGLU>
    F.iref = iref;
end
F.a = a./a(iref);
F.lref  = F.l(1:m)./F.l(iref);
F.lrefc = cumsum(F.lref);
F.C0    = F.C0(1:m);
F.l     = F.l(1:m);
F.D     = F.D(1:m);

% check if restrictions applied %added 24/01/07
if isfield(F,'ivalid')
    F.ivalid = F.ivalid(F.ivalid<m);
    F.a = F.a(F.ivalid);
    F.lref = F.lref(F.ivalid);
    F.lrefc = F.lrefc(F.ivalid);
    F.C0 = F.C0(F.ivalid);
    F.l = F.l(F.ivalid);
    F.D = F.D(F.ivalid);
    F.k = F.k(F.ivalid);
    m = length(F.ivalid);
    F.m = m;
end

% equilibrium value
C0eq = sum(F.lref*F.L.*F.C0)/(1+sum((F.k0./F.k .* F.lref*F.L)));
F.peq = F.k0 * C0eq;

% mesh generation
X = ones(F.m,1);
for i=2:F.m
    X(i) = X(i-1)*(F.a(i-1)*F.lref(i))/(F.a(i)*F.lref(i-1));
end
X = max(F.nmeshmin,ceil(F.nmesh*X/sum(X)));
X = round((X/sum(X))*F.nmesh); % update
F.nmesh = sum(X); % roundoff
xmesh    = zeros(F.nmesh,1);
D        = xmesh;
k        = xmesh;
de       = xmesh; % distance to the next interface in the east direction
dw       = xmesh; % distance to the next interface in the west direction
C0       = xmesh;
j = 1; x0 = 0;
for i=1:F.m
    ind = j+(0:X(i)-1);
    de(ind) = F.lref(i)/(2*X(i));
    dw(ind) = F.lref(i)/(2*X(i));
    D(ind)  = F.D(i)/F.D(F.iref);
    k(ind)  = F.k(i);
    C0(ind) = F.C0(i)/C0eq;
    xmesh(ind) = linspace(x0+dw(ind(1)),F.lrefc(i)-de(ind(end)),X(i));
    x0 = F.lrefc(i);
    j = ind(end)+1;
end

% use a previous solution (if any)
if isfield(F,'restart') && ~isempty(F.restart) && isstruct(F.restart) && isfield(F.restart,'x') && isfield(F.restart,'C')
    if ~isfield(F.restart,'method'), F.restart.method = 'linear'; end
    C0 = interp1(F.restart.x/F.restart.x(end),F.restart.C,xmesh/xmesh(end),F.restart.method)/C0eq; 
end

% equivalent conductances
he = zeros(F.nmesh,1); hw = he;
hw(1) = 1/( (F.k0/k(1))/(F.Bi) + dw(1)/(D(1))   ); % pervious contact
for i=2:F.nmesh
    hw(i) = 1/( (de(i-1)/D(i-1))*(k(i-1)/k(i)) + dw(i)/D(i)  );
end
% for i=1:F.nmesh-1
%     he(i) = hw(i+1); %flux continuity
% end
he(1:F.nmesh-1) = hw(2:F.nmesh); %flux continuity
he(end) = 0; % impervious BC
% control (for debug)
% clf, hold on
% plot(xmesh,zeros(size(xmesh)),'ro')
% plot(xmesh+de,zeros(size(xmesh)),'^')
% plot(xmesh-dw,zeros(size(xmesh)),'v')
% line(repmat(F.lrefc,2,1),repmat(ylim',1,F.m),'color','k')

% Assembling
A = zeros(F.nmesh+1,F.nmesh+1);
% fluid phase = node 0 (position 1)
A(1,1:2) = F.L*hw(1) * [-F.k0/k(1)
                         1           ]';
% node 1 (position 2)
A(2,1:3) = 1/(dw(1)+de(1)) * [ hw(1)*F.k0/k(1)
                                      -hw(1)-he(1)*k(1)/k(2)
                                       he(1) ]';                     
% node i<n (position i+1)
for i=2:F.nmesh-1
    A(i+1,i:i+2) = 1/(dw(i)+de(i)) * [ hw(i)*k(i-1)/k(i)
                                      -hw(i)-he(i)*k(i)/k(i+1)
                                       he(i) ]';
end
% node n (position n+1)
i = F.nmesh;
A(end,end-1:end) = 1/(dw(i)+de(i)) * [ hw(i)*k(i-1)/k(i)
                                      -hw(i)]';
                                                               
% integration
dCdt = @(t,C)(A*C);

function rate = dCdt_opt(t,C) %#ok<INUSL>
rate  = [ ...
    A{1}(1)*C(1) + A{2}(1)*C(2) ;...
    A{1}(2:end-1).*C(2:end-1) + A{2}(2:end).*C(3:end) + A{3}(1:end-1).*C(1:end-2) ;...
    A{1}(end)*C(end) + A{3}(end)*C(end-1) ...
    ];
end

% J = @(t,C)(A);
% F.options.Jacobian = J;
if F.nmesh<400
    F.options.Vectorized = 'on';
    [t,C] = ode15s(dCdt,F.t,[0;C0],F.options); % integration
else % accelerated procedure
    F.options.Vectorized = 'off';
    A = {diag(A,0) diag(A,1) diag(A,-1)};
    [t,C] = ode15s(@dCdt_opt,F.t,[0;C0],F.options); % integration
end

CF = C(:,1);
C  = C(:,2:end);

% interpolation at each interface
Ce = zeros(length(t),F.nmesh);
Cw = zeros(length(t),F.nmesh);
for i=1:length(t)
    Ce(i,1:end-1) = C(i,1:end-1) - (de(1:end-1).*he(1:end-1).*( k(1:end-1)./k(2:end) .* C(i,1:end-1)' - C(i,2:end)' ) ./ D(1:end-1) )';
    Ce(i,end)     = C(i,end-1);
    Cw(i,2:end)   = C(i,2:end)   + (dw(2:end)  .*hw(2:end)  .*( k(1:end-1)./k(2:end) .* C(i,1:end-1)' - C(i,2:end)' ) ./ D(2:end)   )';
    Cw(i,1)       = C(i,1)       +  dw(1)       *hw(1)       *( F.k0       /k(1)      * CF(i)         - C(i,1)      )  / D(1);
end
Cfull = reshape([Cw;C;Ce],length(t),3*F.nmesh);
xw = xmesh-dw+zero;
xe = xmesh+de-zero;
xfull = [xw';xmesh';xe']; xfull = xfull(:);
kfull = [k';k';k']; kfull = kfull(:); %#ok<NASGU>

% outputs
res.C = interp1(t,trapz(xfull,Cfull,2)*C0eq,ti,method)/xfull(end); % av conc
res.t = ti; % time base
res.p = NaN; %interp1(t,repmat(kfull',length(t),1).*Cfull*C0eq,ti,method); % to fast calculations
res.peq = F.peq;
res.V  = F.lrefc(end); %*F.l(F.iref);
res.C0 = C0*C0eq;
res.F  = F;
res.x  = xfull; %*F.l(F.iref);
res.Cx = Cfull*C0eq; %interp1(t,Cfull*C0eq,ti,method);
res.tC = t;
res.CF = interp1(t,CF*C0eq,ti,method);
res.fc = res.CF/F.L;
res.f  = interp1(t,hw(1) * ( F.k0/k(1) * CF - C(:,1) ) * C0eq,ti,method);
res.timebase = res.F.l(res.F.iref)^2/(res.F.D(res.F.iref));

end

% example
% F	= 	struct(...
% 				'Bi'		, 	1e3,...	Biot [hm.L1/D]
% 				'k'			,	[5 1],...[0.5 3 2],...	ki, i=1 (layer in contact with the liquid)
%                 'D'         ,   [1e-14 1e-14],... diffusion coefficient
%                 'k0'        ,   1,... 0 = liquid
%                 'l'         ,   [50 100]*1e-6,...[50 20 10 120]*1e-6,... m
% 				'L'			,	1/10,...	dilution factor (respectively to iref)
% 				'C0'		,	[0 500],...	initial concentration in each layer
%                 't'         ,   [0 10] ...
% 					); % if iref is mission, it is indentified
% day=24*3600;
% tunit = res.timebase/day;
% tbase = res.t*tunit;
% tsamp = res.tC*tunit;
% plot(tbase,res.CF)
% plot(res.x,res.Cx')