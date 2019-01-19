function hout = plotsmoothbar(x,y,removeflat,varargin)
%  PLOTSMOOTHBAR plots a smoothed bar centered on X
%   h=plotsmoothbar(x,y,removeflat,[...])
%       [...] = pair properties as plot

% MS-MATLAB 1.0 - 28/06/07 - Olivier Vitrac - rev. 26/04/2018

% revision history
% 26/04/2018 compatibility with R2014b and later

% definitions
method = 'pchip';
amplification = 30;
removeflat_default = 0;

% arg check
if nargin<2, error('2 arguments are at least required'), end
if nargin<3, removeflat = []; end
if isempty(removeflat), removeflat = removeflat_default; end

% size
dx = max(diff(x));
m  = length(x);
n  = amplification*m+1;
c  = size(y,2);
xi = linspace(x(1)-dx/2,x(end)+dx/2,n)';

% plot
yi = interp1(x,y,xi,method);
if removeflat
    nz = yi>sqrt(eps);
    ind = nz(min(n,find(abs(diff(yi(nz)))<sqrt(eps))+1));
    nind = setdiff(1:n,ind);
    yi(ind) = interp1(xi(nind),yi(nind),xi(ind),method); %,'extrap');
end
yi(yi<0)=0;
yi = yi./repmat(trapz(xi,yi,1),n,1); %repmat(sum(y,1),n,1) .* 
h = plot(xi,yi,varargin{:});

if nargout, hout = h; end
