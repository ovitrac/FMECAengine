function [xx,yy,zz] = cylrevolution(r,n,sq,rmax,heights,xscaling,yscaling)
%CYLREVOLUTION Generate a 3D object by cylindrical revolution cylinder.
%   [X,Y,Z] = CYLREVOLUTION(R[,N,ISQUARE,RMAX,xscaling,yscaling,zscaling]) forms the unit cylinder based on the generator
%   curve in the vector R. Vector R contains the radius at equally
%   spaced points along the unit height of the cylinder. The cylinder
%   has N points around the circumference. SURF(X,Y,Z) displays the
%   cylinder.
%
%   Omitting output arguments causes the cylinder to be displayed with
%   a SURF command and no outputs to be returned.
%
%   CYLREVOLUTION(AX,...) plots into AX instead of GCA.
%
%   See also CYLINDER, SPHERE, ELLIPSOID.
%
%   Example:
%       cylrevolution([1 2 0.5],300,[true,false,true],[1.3 Inf])
%{
   % Advanced example
     cylrevolution( [0, 0.1, 0.3, 1, 1 0.25 0.25 0.22 0.22],...
                    300,...
                    [false, false, false, true, true, false, false false false],...
                    [Inf, Inf  , Inf  , 1.3 , 1.3, Inf, Inf Inf Inf],...
                    [0, 0 , -0.05 , 0   , 3, 0.05, 0.3, 0 -0.33], ...
                    [1, 1 , 1     , 0.5   , 0.5, 1   , 1  , 1 1] ...
                  );
%}
%{
   % Advanced example
   B = array2table([
         0         0       Inf         0    1.0000   % - conical bottom (outer)
    0.1000         0       Inf         0    1.0000
    0.3000         0       Inf   -0.0500    1.0000   
    1.0000    1.0000    1.3000         0    0.5000   % - vetical walls
    1.0000    1.0000    1.3000    3.0000    0.5000
    0.2500         0       Inf    0.0500    1.0000
    0.2500         0       Inf    0.3000    1.0000
    0.2200         0       Inf         0    1.0000
    0.2200         0       Inf   -0.3300    1.0000
    0.97           1      1.27   -0.05      0.5
    0.97           1      1.27   -2.97     0.5
    0.3            0      Inf    0.053      1
    0.1            0      Inf    0          1
    0              0      Inf    0          1
    ],'VariableNames',{'r' 'isquare' 'rmax' 'H' 'W'});
    cylrevolution(B.r,300,B.isquare,B.rmax,B.H,B.W)
%}
%
% based on Matlab function CYLINDER

% MS 2.1 - 16/12/2017 - INRA\Olivier Vitrac - rev. 

% default
r_default = [1;1];
n_default = 20;
rmax_default = Inf;
heights_default = 1;
xscaling_default = 1;
yscaling_default = 1;

% arg check
if nargin<1, r = []; end
if nargin<2, n = []; end
if nargin<3, sq = []; end
if nargin<4, rmax = []; end
if nargin<5, heights = []; end
if nargin<6, xscaling = []; end
if nargin<7, yscaling = []; end
if isempty(r), r = r_default; end, r = r(:);
if isempty(n), n = n_default; end
[cax,args,nargs] = axescheck(r,n);
if nargs > 0, r = args{1}; end
if nargs > 1, n = args{2}; end
if isempty(sq), sq = true(size(r)); end, sq = sq>0; sq = sq(:);
if isempty(rmax), rmax = rmax_default; end, rmax = rmax(:);
if isempty(heights), heights = heights_default; end, heights = heights(:);
if isempty(xscaling), xscaling = xscaling_default; end, xscaling = xscaling(:);
if isempty(yscaling), yscaling = yscaling_default; end, yscaling = yscaling(:);


% preparation
r(r==0) = eps;
m = length(r); if m==1, r = [r;r]; m = 2; end
sq = argpad(sq,m+1); rmax = argpad(rmax,m); sq = argpad(sq,m);
theta = (0:n)/n*2*pi;
sintheta = sin(theta); sintheta(n+1) = 0;
sqfunc = ones(m,1) * min(1./abs(cos(theta)),1./abs(sintheta)); sqfunc(~sq,:) = 1;
rmax = rmax*ones(size(theta));
heights = argpad(heights,m); heights = cumsum(heights)-heights(1);
xscaling = argpad(xscaling,m); yscaling = argpad(yscaling,m); 

% plots
O = ones(1,n+1);
x = (r * cos(theta)) .*sqfunc;
y = (r * sintheta  ) .*sqfunc;
rcalc = sqrt(x.^2+y.^2);
rscale = min(rmax,rcalc)./rcalc;
x = x.*rscale.*(xscaling*O);
y = y.*rscale.*(yscaling*O);
z = heights * O;

if nargout == 0
    cax = newplot(cax);
    surf(x,y,z,'parent',cax,'edgecolor','none','facealpha',.5)
    camlight left, lighting gouraud, axis equal
else
    xx = x; yy = y; zz = z;
end
