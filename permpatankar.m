function res = permpatankar(varargin)
%SENSPATANKAR simulates transfer through n layers using a modified Patankar Method (see p 45)
%   Input constructor :
%       S = permpatankar()
%       S = permpatankar('init' [,'property1',value1,'property2',value2,...])
%       S = permpatankar('default' [,'property1',value1,'property2',value2,...])
%   Run simulation S
%       R = permpatankar(S);
%   Supersede parameters in S before launching simulation
%       R = permpatankar('property1',value1,'property2',value2,...,S)
%
%   List of implemented properties/values (default, here 3 layers)
%                  Bi: mass Biot number in the upstream/downstream ([1000 1000])
%                   k: Henry-like coefficients for the three layers ([1 1 1])
%                   D: diffusion coefficients ([1.0000e-14 1.0000e-14 1.0000e-14])
%                  k0: Henry-like coefficients for the upstream and downstream compartments [1 1]
%                   l: layer thickness [3.0000e-05 3.0000e-05 3.0000e-05]
%                   L: dilution ratios lupstream/sum(li) and ldownstream/sum(li) ([0.0100 0.1000])
%                  C0: initial concentration within each layer ([500 500 500])
%                 CF0: initial concentration in the upstream and downstream compartment ([0 0])
%                iref: reference layer for dimension formulation ([])
%              ploton: plot flag (false)
%              dispon: disp flag (false)
%               debug: debug flag (false)
%            autotime: 'on'
%                   n: 10000
%               nmesh: number of internal nodes (200)
%            nmeshmin: number of minimum nodes by layer (20)
%              method: interpolation method ('pchip')
%          maxstepmax: maximum integration time steps (100)
%         nstepchoice: depreciated, not used (200)
%                zero: zero value (2.2204e-14)
%                   t: Fourier time when to get the solution ([1x560 double], [0:.00001:.005 .01:.01:.1 .11:.1:5])
%
%       OUTPUT with fields
%                C: [10000x1 double]
%                t: [10000x1 double]
%                p: NaN
%              peq: 909.0884
%                V: 1
%               C0: [200x1 double]
%                F: [1x1 struct]
%                x: [600x1 double]
%               Cx: [560x600 double]
%               tC: [560x1 double]
%              CF1: [10000x1 double]
%              CF2: [10000x1 double]
%            CF1eq: 909.0884
%            CF2eq: 909.0884
%              fc1: [10000x1 double]
%              fc2: [10000x1 double]
%               f1: [10000x1 double]
%               f2: [10000x1 double]
%         timebase: 9.0000e+08
%       
%
% Basic example based on default properties
%{
  S = permpatankar();
  R = permpatankar(S);
%}
%
% Example with over-definitions
%{
  S0 = permpatankar('init');
  R=permpatankar('k',1,'D',1e-16,'l',300e-6,'C0',0,'CF0',[1000 0],S0);
  figure, subplot(2,1,1), plot(R.t,R.CF1), subplot(2,1,2), plot(R.t,R.CF2)
%}


% MS-MATLAB-WEB 1.0 - 28/09/07 - Olivier Vitrac - rev 29/01/2016

% Revision history
% 25/12/2016 better compatibility with R2016, new arguments control (to be further checked before being distributed)
% 03/01/2017 fix LATeX equations, fix confusion between l and L, improved help, new constructor
% 26/01/2017 separation of l and L has been removed
% 29/01/2017 fix cumulative flux

% definitions
global timeout
timeout		= 800; % s
defaultoptions  = odeset('RelTol',1e-4,'AbsTol',1e-4,'Initialstep',1e-8,'Maxstep',.01,'Maxorder',2); %odeset('RelTol',1e-4,'AbsTol',1e-4,'Stats','yes','Initialstep',1e-5,'Maxstep',.05,'Maxorder',5);
Fdefault	= 	struct(...
				'Bi'		, 	[1e3 1e3],...	Biot [hm.L1/D]
				'k'			,	[1 1 1],...[0.5 3 2],...	ki, i=1 (layer in contact with the liquid)
                'D'         ,   [1e-14 1e-14 1e-14],... diffusion coefficient
                'k0'        ,   [1 1],... 0 = liquid
                'l'         ,   [30 30 30]*1e-6,...[50 20 10 120]*1e-6,... m
				'L'			,	[.01 .1],...	dilution factor (respectively to iref)
				'C0'		,	[500 500 500 ],...	initial concentration in each layer
                'CF0'       ,   [0 0],...
                'iref'      , [],...
                'ploton',false,...
                'dispon',false,...
                'debug',false,...
                'autotime','on',...
                'n',1e4,...
                'nmesh',200,... % number of nodes for a layer of normalized thickness 1
                'nmeshmin',20,...
                'method','pchip',...
                'maxstepmax',100,...
                'nstepchoice',200,...
                'zero',100*eps ...    
					); % if iref is missing, it is indentified
                
% Default dimensionless times to calculate the solution
Fdefault.t	= [0:.00001:.005 .01:.01:.1 .11:.1:5]; %0:.00001:.005; %[0:.00001:.005 .01:.01:.1 .11:.1:5]';
if min(Fdefault.Bi)<10, Fdefault.t = Fdefault.t/(min(Fdefault.Bi)/10); end

% arg check
if ~nargin
    initon = true;
elseif ischar(varargin{1})
    if isempty(regexpi(varargin{1},'^default|^init.*'))
       initon = false;
    else
        initon = true;
        if nargin>1
            Fraw = argcheck(varargin{2:end},Fdefault,'','keep','nostructexpand','case');
            if isfield(Fraw,'options')
            options = argcheck(Fraw.options,defaultoptions,'');
            else
                options = defaultoptions;
            end
            Fdefault = argcheck(varargin{2:end},Fdefault,'','case');
            Fdefault.options = options;
        end
    end
else
    initon = false;
end
if ~initon % extract options
    Fraw = argcheck(varargin,Fdefault,'','keep','nostructexpand','case'); % keep case, structures and extra parameters
    if isfield(Fraw,'options')
        options = argcheck(Fraw.options,defaultoptions,'');
    else
        options = defaultoptions;
    end
    F = argcheck(varargin,Fdefault,'','case'); % use 'case' rather than a separation of l and L
    % restore the separation between l and L
    %F.l = Fraw.l;
    %F.L = Fraw.L;
else
    options = defaultoptions;
    F = Fdefault;
end


%% PHYSICAL CHECK
% check that properties are positive
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

% returns default if initizalization is true
if initon
    res = Fdefault;
    if F.dispon, disp('... init mode'), end
    return
end

if strcmpi(F.autotime,'on')
	ti		= linspace(min(F.t),max(F.t),F.n)';
else
	ti		= F.t;
end
F.maxstepmax = max(F.maxstepmax,ceil(F.t(end)/F.nstepchoice));

%% INPUT RENORMALIZATION (fields X1->Xn)
% iref is the index of the reference barrier (with the least permeability D/(k*l)
F.k = F.k/F.k0(1); F.k0 = F.k0/F.k0(1);   % by convention
a = F.D(1:m)./F.k(1:m);
if isfield(F,'iref') && ~isempty(F.iref)
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

%% EQUILIBRIUM CONCENTRATION in upstream compartment
% 
% $${C_{{F_{eq}},1}} = \frac{{\frac{{{C_{{F_0},1}}}}{{{L_1}}} + \frac{{{C_{{F_0},2}}}}{{{L_2}}} + \sum\limits_{i = 1}^{{n_{layer}}} {\frac{{{l_i}}}{{{l_{{i_{ref}}}}}}{C_{0,i}}} }}{{\frac{1}{{{L_1}}} + \frac{{{k_{0,1}}}}{{{k_{0,2}}}}\frac{1}{{{L_2}}} + \sum\limits_{i = 1}^{{n_{layer}}} {\frac{{{k_{0,1}}}}{{{k_i}}}\frac{{{l_i}}}{{{l_{{i_{ref}}}}}}} }}$$
% 
massbalance = F.CF0(1)/F.L(1) + F.CF0(2)/F.L(2) + sum( (F.l/F.lref).*F.C0 );
eq1 =  1/F.L(1) + (F.k0(1)/F.k0(2))*1/F.L(2) + sum(  (F.k0(1)./F.k) .* (F.l/F.lref) ); 
CFeq1 = massbalance / eq1;
%% EQUILIBRIUM CONCENTRATION in downstream compartment
% 
% $${C_{{F_{eq}},2}} = \frac{{\frac{{{C_{{F_0},1}}}}{{{L_1}}} + \frac{{{C_{{F_0},2}}}}{{{L_2}}} + \sum\limits_{i = 1}^{{n_{layer}}} {\frac{{{l_i}}}{{{l_{{i_{ref}}}}}}{C_{0,i}}} }}{{\frac{1}{{{L_2}}} + \frac{{{k_{0,2}}}}{{{k_{0,1}}}}\frac{1}{{{L_1}}} + \sum\limits_{i = 1}^{{n_{layer}}} {\frac{{{k_{0,2}}}}{{{k_i}}}\frac{{{l_i}}}{{{l_{{i_{ref}}}}}}} }}$$
% 
eq2 =  1/F.L(2) + (F.k0(2)/F.k0(1))*1/F.L(1) + sum(  (F.k0(2)./F.k) .* (F.l/F.lref) ); 
CFeq2 = massbalance / eq2;
%% NORMALIZATION VALUE
C0eq = ( CFeq1/F.L(1)+CFeq2/F.L(2) ) / sum(1./F.L);
% C0eq = (    sum( F.lref*F.L(1).*F.C0) + ...
%             1*F.CF0(1) + ...
%             F.L(2)/F.L(1)*F.CF0(2) )/ ...
%         (1 + F.k0(1)/F.k0(2)*F.L(2)/F.L(1) + sum((F.k0(1)./F.k .* F.lref*F.L(1))));
F.peq = F.k0(1) * CFeq1;

%% MESH GENERATION
% Create 1D mesh (volume interfaces and internal nodes), assign properties (uniform) within each volume
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

if F.debug
    figure, hold on
    hp = [plot(xmesh,zeros(size(xmesh)),'ro')
          plot(xmesh+de,zeros(size(xmesh)),'^')
          plot(xmesh-dw,zeros(size(xmesh)),'v')
         line(repmat(F.lrefc,2,1),repmat(ylim',1,F.m),'color','k')
         ];
    hl=legend(hp(1:4),{'\bfnode' 'w' 'e'},'fontsize',10,'location','best'); set(hl,'box','off')
    formatax(gca,'fontsize',10)
    title(sprintf('%d layers',F.m),'fontsize',12)
end

%% EQUIVALENT CONDUCTANCES
% conductances are defined at the interface with the upstream compartment
he = zeros(F.nmesh,1); hw = he;

%% conductance at the interface with the upstream compartment
% 
% $${h_{w,1}} = \frac{1}{{\frac{{{k_{0,1}}}}{{{k_1}B{i_1}}} + \frac{{{d_{w,1}}}}{{{D_1}}}}}$$
% 
hw(1) = 1/( (F.k0(1)/k(1))/(F.Bi(1)) + dw(1)/(D(1))   ); % pervious contact
%% conductance for internal nodes
% 
% $${h_{w,i}} = \frac{1}{{\frac{{{d_{e,i - 1}}}}{{{D_{i - 1}}}}\frac{{{k_{i - 1}}}}{{{k_i}}} + \frac{{{d_{w,i}}}}{{{D_i}}}}}$$
% 
for i=2:F.nmesh
    hw(i) = 1/( (de(i-1)/D(i-1))*(k(i-1)/k(i)) + dw(i)/D(i)  );
end
%% continuity equation
% 
% $${h_{e,1..n - 1}} = {h_{w,2..n}}$$
% 
he(1:F.nmesh-1) = hw(2:F.nmesh); %flux continuity
%% conductance at the interface with the upstream compartment
% 
% $${h_{e,n}} = \frac{1}{{\frac{1}{{B{i_2}}} + \frac{{{k_n}}}{{{k_{0,2}}}}\frac{{{d_{e,n}}}}{{{D_n}}}}}$$
% 
he(end) = 1/( 1/F.Bi(2) + (k(end)/F.k0(2))*de(end)/(D(end))   ); % pervious BC


%% RIGIDITY MATRIX
% (n+2)x(n+2) matrix codes for the mass balance. First and last row code for the mass balance in the upstream and downstreal compartment
A = zeros(F.nmesh+2,F.nmesh+2);

%% Upstream fluid node (position 1)
% 
% $$\frac{{d{C_1}}}{{dt}} = {L_1}{h_{w,1}}\left( { - \frac{{{k_{0,1}}}}{{{k_1}}}{C_1} + {C_2}} \right) =  = {L_1}{h_{w,1}}\left( {{C_2} - \frac{{{k_{0,1}}}}{{{k_1}}}{C_1}} \right) = {L_1}\frac{{{C_2} - \frac{{{k_{0,1}}}}{{{k_1}}}{C_1}}}{{\frac{{{k_{0,1}}}}{{{k_1}B{i_1}}} + \frac{{{d_{w,1}}}}{{{D_1}}}}} = {L_1}\frac{{\frac{{{k_1}}}{{{k_{0,1}}}}{C_2} - {C_1}}}{{\frac{1}{{B{i_1}}} + \frac{{{k_1}}}{{{k_{0,1}}}}\frac{{{d_{w,1}}}}{{{D_1}}}}}$$
% 
A(1,1:2) = F.L(1)*hw(1) * [-F.k0(1)/k(1)
                         1           ]';
%% node 1 (position 2)
%
% $$\frac{{d{C_2}}}{{dt}} = \frac{1}{{{d_{w,1}} + {d_{e,1}}}}\left( {{h_{w,1}}\frac{{{k_{0,1}}}}{{{k_1}}}{C_1} - \left( {{h_{w,1}} + {h_{e,1}}\frac{{{k_1}}}{{{k_2}}}} \right){C_2} + {h_{e,1}}{C_3}} \right)$$
%
A(2,1:3) = 1/(dw(1)+de(1)) * [ hw(1)*F.k0(1)/k(1)
                                      -hw(1)-he(1)*k(1)/k(2)
                                       he(1) ]';                     
%% internal nodes i<n (position i+1)
% 
% $$\frac{{d{C_{i + 1}}}}{{dt}} = \frac{1}{{{d_{w,i}} + {d_{e,i}}}}\left( {{h_{w,i}}\frac{{{k_{i - 1}}}}{{{k_i}}}{C_i} - \left( {{h_{w,i}} + {h_{e,i}}\frac{{{k_i}}}{{{k_{i + 1}}}}} \right){C_{i + 1}} + {h_{e,i}}{C_{i + 1}}} \right)$$
% 
for i=2:F.nmesh-1
    A(i+1,i:i+2) = 1/(dw(i)+de(i)) * [ hw(i)*k(i-1)/k(i)
                                      -hw(i)-he(i)*k(i)/k(i+1)
                                       he(i) ]';
end

%% terminal node n (position n+1)
% 
% $$\frac{{d{C_{n + 1}}}}{{dt}} = \frac{1}{{{d_{w,n}} + {d_{e,n}}}}\left( {{h_{w,n}}\frac{{{k_{n - 1}}}}{{{k_n}}}{C_n} - \left( {{h_{w,n}} + {h_{e,n}}\frac{{{k_n}}}{{{k_{0,2}}}}} \right){C_{n + 1}} + {h_{e,n}}{C_{n + 2}}} \right)$$
% 
i = F.nmesh;
A(i+1,i:i+2) = 1/(dw(i)+de(i)) * [ hw(i)*k(i-1)/k(i)
                                      -hw(i)-he(i)*k(i)/F.k0(2)
                                       he(i) ]';
                                   
%% fluid phase (node n+1, position n+2)
% 
% $$\frac{{d{C_{n + 2}}}}{{dt}} = {L_2}{h_{e,n}}\left( {\frac{{{k_{_{n + 2}}}}}{{{k_{0,2}}}}{C_{n + 1}} - {C_{n + 2}}} \right) = {L_2}\frac{{\frac{{{k_{_{n + 2}}}}}{{{k_{0,2}}}}{C_{n + 1}} - {C_{n + 2}}}}{{\frac{1}{{B{i_2}}} + \frac{{{k_n}}}{{{k_{0,2}}}}\frac{{{d_{e,n}}}}{{{D_n}}}}}$$
% 
A(end,end-1:end) = F.L(2)*he(end) * [F.k(end)/F.k0(2)
                                     -1 ]';
                                 
%% INTEGRATION OF TRANSPORT EQUATION
% linearized and vectorized version: dCdt = @(t,C)(A*C);
dCdt = @(t,C)(A*C);
% J = @(t,C)(A);
% F.options.Jacobian = J;
options.Vectorized = 'on'; % 'on' only if vectorized
[t,C] = ode15s(dCdt,F.t,[F.CF0(1)/C0eq;C0;F.CF0(2)/C0eq],options); % integration
CF1 = C(:,1);
CF2 = C(:,end);
C  = C(:,2:end-1);

% interpolation at each interface
Ce = zeros(length(t),F.nmesh);
Cw = zeros(length(t),F.nmesh);
for i=1:length(t)
    Ce(i,1:end-1) = C(i,1:end-1) - (de(1:end-1).*he(1:end-1).*( k(1:end-1)./k(2:end) .* C(i,1:end-1)' - C(i,2:end)' ) ./ D(1:end-1) )';
    Ce(i,end) =     C(i,end) - (de(end) * he(end).*( k(end)./F.k0(2) .* C(i,end) - CF2(i) ) / D(end) );
    Cw(i,2:end)   = C(i,2:end)   + (dw(2:end)  .*hw(2:end)  .*( k(1:end-1)./k(2:end) .* C(i,1:end-1)' - C(i,2:end)' ) ./ D(2:end)   )';
    Cw(i,1)       = C(i,1)       +  dw(1)       *hw(1)       *( F.k0(1)       /k(1)      * CF1(i)         - C(i,1)      )  / D(1);
end
Cfull = reshape([Cw;C;Ce],length(t),3*F.nmesh);
xw = xmesh-dw+F.zero;
xe = xmesh+de-F.zero;
xfull = [xw';xmesh';xe']; xfull = xfull(:);
kfull = [k';k';k']; kfull = kfull(:); %#ok<NASGU>

% outputs
res.C = interp1(t,trapz(xfull,Cfull,2)*C0eq,ti,F.method)/xfull(end); % av conc
res.t = ti; % time base
res.p = NaN;
res.peq = F.peq;
res.V  = F.lrefc(end); %*F.l(F.iref);
res.C0 = C0*C0eq;
res.F  = F;
res.x  = xfull; %*F.l(F.iref);
res.Cx = Cfull*C0eq; %interp1(t,Cfull*C0eq,ti,method);
res.tC = t;
res.f1  = interp1(t,hw(1) * ( F.k0(1)/k(1) * CF1 - C(:,1) ) * C0eq,ti,F.method);
res.f2  = interp1(t,he(end) * ( F.k(end)/F.k0(2) * C(:,end) - CF2 ) * C0eq,ti,F.method);
res.CF1 = interp1(t,CF1*C0eq,ti,F.method);
res.CF2 = interp1(t,CF2*C0eq,ti,F.method);
res.fc1 = (res.CF1-res.CF1(1))*F.l(F.iref)/F.L(1);
res.fc2 = (res.CF2-res.CF2(1))*F.l(F.iref)/F.L(2);
res.CF1eq = CFeq1;
res.CF2eq = CFeq2;

res.timebase = res.F.l(res.F.iref)^2/(res.F.D(res.F.iref));

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