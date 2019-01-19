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

% MS 2.1 - 16/12/2017 - INRA\Olivier Vitrac - rev. 28/10/2018

% Revision history
% 27/10/2018 implement differences between left and right
% 28/10/2018 add r as a structure with fields, r, erodefunc() and smooth
% 29/10/2018 RC with r as a structure

% default
r_default = [1;1];
n_default = 20;
rmax_default = Inf;
heights_default = 1;
xscaling_default = 1;
yscaling_default = 1;
erodefunc = [];
smooth = 0;

% arg check
if nargin<1, r = []; end
if nargin<2, n = []; end
if nargin<3, sq = []; end
if nargin<4, rmax = []; end
if nargin<5, heights = []; end
if nargin<6, xscaling = []; end
if nargin<7, yscaling = []; end
if isstruct(r)
    if ~isfield(r,'erodefunc'), error('the field erodefun is missing in r object'), end
    if ~isfield(r,'erodefunc'), error('the field r is missing in r object'), end
    erodefunc = r.erodefunc;
    if ~isa(erodefunc,'function_handle') || isempty(erodefunc), error('erodefun should be empty or a function_handle'), end
    if isfield(r,'smooth'), smooth = r.smooth; end
    if isfield(r,'r'), r = r.r; end
    if numel(smooth)>1, error('smooth must be a scalar value'), end
    if smooth>=1, error('smooth (current value = %0.3g) cannot be larger than unity',smooth), end
end
if isempty(r), r = r_default; end, r = r(:);
if isempty(n), n = n_default; end
[cax,args,nargs] = axescheck(r,n);
if nargs > 0, r = args{1}; end
if nargs > 1, n = args{2}; end
if isempty(sq), sq = true(size(r)); end, sq = sq>0; if ndims(sq)<2, sq = sq(:); end
if isempty(rmax), rmax = rmax_default; end, if ndims(rmax)<2, rmax = rmax(:); end
if isempty(heights), heights = heights_default; end, heights = heights(:);
if isempty(xscaling), xscaling = xscaling_default; end, xscaling = xscaling(:);
if isempty(yscaling), yscaling = yscaling_default; end, yscaling = yscaling(:);

% preparation
r(r==0) = eps;
m = length(r); if m==1, r = [r;r]; m = 2; end
%sq = argpad(sq,m+1);
itmp = argpad(1:size(sq,1),m); sq =sq(itmp,:);
itmp = argpad(1:size(rmax,1),m); rmax = rmax(itmp,:);
theta = (0:n)/n*2*pi;sintheta = sin(theta); sintheta(n+1) = 0;
sqfunc = ones(m,1) * min(1./abs(cos(theta)),1./abs(sintheta));
isatleft = ((ones(m,1) * cos(theta))<0);
if size(sq,2)==1
    sqfunc(~sq,:) = 1;
else
    [sqfuncleft,sqfuncright] = deal(sqfunc);
    sqfuncleft(~sq(:,1),:) = 1;
    sqfuncright(~sq(:,2),:) = 1;
    sqfunc(isatleft) = sqfuncleft(isatleft);
    sqfunc(~isatleft) = sqfuncright(~isatleft);    
end
heights = argpad(heights,m); heights = cumsum(heights)-heights(1);
xscaling = argpad(xscaling,m); yscaling = argpad(yscaling,m); 

% plots
O = ones(1,n+1);
x = (r * cos(theta)) .*sqfunc;
y = (r * sintheta  ) .*sqfunc;
rcalc = sqrt(x.^2+y.^2);
if size(rmax,2)==1
    rmax = rmax*ones(size(theta));
    rscale = min(rmax,rcalc)./rcalc;
    ratiox = 1;
else
    ratiox = (max(rmax(~isinf(rmax(:,2)),2)) / max(rmax(~isinf(rmax(:,1)),1)));
    rcalc(~isatleft) = rcalc(~isatleft);
    rmax1 = rmax(:,1)*ones(size(theta));
    rmax2 = rmax(:,2)*ones(size(theta));
    rscale = rmax1;
    rscale(isatleft) = min(rmax1(isatleft),rcalc(isatleft))./rcalc(isatleft);
    rscale(~isatleft) = min(rmax2(~isatleft),rcalc(~isatleft))./rcalc(~isatleft);
end
x = x.*rscale.*(xscaling*O);
if (size(rmax,2)>1) && (size(sq,2)>1)
    irawtofix = find(~isinf(rmax(:,2)) & (sq(:,2)==1));
    xtmp = x(irawtofix,:);
    ytmp = y(irawtofix,:);
    if ~isempty(erodefunc)
        xmax = max(xtmp(:));
%         anglefunc = @(scale,x,y) (y(:)-0) ./ (x(:)-2.8*scale);
%         cosfunc = @(scale,x,y) scale*(2.8-2*abs(cos(anglefunc(scale,x,y))));
%         sinfunc = @(scale,x,y) 2*scale*sign(y(:)).*abs(sin(anglefunc(scale,x,y)));
%         erodefunc = @(scale,x,y) [cosfunc(scale,x,y) sinfunc(scale,x,y)];
        xy = erodefunc(xmax,xtmp,ytmp);
        xtmpi = interpleft(xy(:,2),xy(:,1),ytmp); % erosion xtmpi, ytmpi
        % envelope
        yscale = linspace(min(ytmp(xtmp>0)),max(ytmp(xtmp>0)),51)';
        xscale = arrayfun(@(i) max(xtmp( (ytmp>=yscale(i)) & (ytmp<yscale(i+1)) )), (1:length(yscale)-1)' );
        xenv = spline(yscale(1:end-1)+diff(yscale)/2,xscale,ytmp);
        % deformation
        dexcess = xtmp - xtmpi; % sqrt(sum(xtmp(:).^2+ytmp(:).^2,2)) - sqrt(sum(xy.^2,2));
        icurve = dexcess>=0;
        translation = dexcess(icurve) + (xenv(icurve)-xtmp(icurve));
        xtmp(icurve) = xtmp(icurve) - translation;
        % smoothing if included
        if smooth>0
            smooth = 0.2;
            nx = size(xtmp,2);
            ibase = unique(round(linspace(1,nx,min(round(nx/10),nx*(1-smooth))))');
            t = linspace(0,1,nx)';
            for jfix = 1:length(irawtofix)
                xytmp= spline(t(ibase),[xtmp(jfix,ibase);ytmp(jfix,ibase)],t);
                xtmp(jfix,:) = xytmp(1,:);
                ytmp(jfix,:) = xytmp(2,:);
            end
        end % endif smooth
    end % endif erodefunc
    xtmp(xtmp>0) = xtmp(xtmp>0) / ratiox;
    
    x(irawtofix,:) = xtmp;
    y(irawtofix,:) = ytmp;
end
y = y.*rscale.*(yscaling*O);
z = heights * O;


if nargout == 0
    cax = newplot(cax);
    surf(x,y,z,'parent',cax,'edgecolor','none','facealpha',.5)
    camlight left, lighting gouraud, axis equal
else
    xx = x; yy = y; zz = z;
end
