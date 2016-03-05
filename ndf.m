function dydt = ndf(t,y,ordre,dydt0,varargin)
%  NDF Numerical differentiation at order 1 with an approximation [h6y(8)/140] or order 2 with an approximation [h6y(8)/560]
%   Syntax: dydt = ndf(t,y [,order,dydt0 [,property,value,...]])
%           t : mx1 array coding for times or x (diff(t) must be constant, no control is done for efficiency) 
%           y : mxn array coding for y values
%       order : 1 (default) or 2
%       dydt0 : initial slope (if known)
%        dydt : first or second derivative in t
%   Pair/property value
%       method: interpolation method (default='pchip')
%  makeuniform: makes data uniform before derivation
%   resolution: size of the uniform grid (default=1e4)

% MS-MATLAB 1.0 - 28/01/04 - Olivier Vitrac - rev. 27/02/16

%Revision history
% 24/05/10 columnwise
% 26/09/15 'cubic' replaced by 'pchip'
% 27/02/16 implement argcheck

% default
default = struct( 'method','pchip','makeuniform',false,'resolution',1e4);
                 
% arg check
if nargin<2, y = t; t = (1:size(y,1))'; end
if nargin<3, ordre = 1; end
if nargin<4, dydt0 = []; end
o = argcheck(varargin,default);
if size(y,1)>1 && size(y,2)>1
    ny = size(y,2);
    dydt = zeros(size(y));
    for i=1:ny
        if any(dydt0)
            dydt(:,i)=ndf(t,y(:,i),ordre,dydt0(min(i,numel(dydt0))),o);
        else
            dydt(:,i)=ndf(t,y(:,i),ordre,[],o);
        end
    end
    return
end

% implements uniform sampling if requested
x0 = [];
if o.makeuniform
    x0 = t;
    if length(unique(x0))==1
        dydt = Inf(size(y));
        return
    end
    t = linspace(min(t),max(t),o.resolution)';
    y = interp1(x0,y,t,o.method);
end
dt	= t(2)-t(1);

switch ordre
case 1
	D	= [-1 +9 -45 0 +45 -9 +1]'/(60*dt);
	D0	= [-137 300 -300 200 -75 12]'/(60*dt);
	D1	= -flipud(D0);
case 2
	D	= [2 -27 +270 -490 270 -27 +2]'/(180*dt.^2);
	D0	= [45 -154 +214 -156 61 -10]'/(12*dt.^2);
	D1 = flipud(D0);
end
% yfull	= [repmat(y(1),3,1) ; y ; repmat(y(end),3,1) ];
% dydt	= [yfull(1:end-6) yfull(2:end-5) yfull(3:end-4) yfull(4:end-3) yfull(5:end-2) yfull(6:end-1)  yfull(7:end)] * D;
dydt = 	[
		[y(1:3) y(2:4) y(3:5) y(4:6) y(5:7) y(6:8) ] * D0;
		[y(1:end-6) y(2:end-5) y(3:end-4) y(4:end-3) y(5:end-2) y(6:end-1)  y(7:end)] * D
		[y(end-7:end-5) y(end-6:end-4) y(end-5:end-3) y(end-4:end-2) y(end-3:end-1)  y(end-2:end)] * D1
	];
if any(dydt0), dydt(1) = dydt0; end
dy = diff(dydt);
i = intersect(find(dy*sign(mean(sign(dy)))<0),2:10);
if any(i), ni = setdiff(1:length(t),i); dydt(i) = interp1(t(ni),dydt(ni),t(i),o.method,'extrap'); end

% back transformation
if o.makeuniform
    dydt = interp1(t,dydt,x0,o.method);
end