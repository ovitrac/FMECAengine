function curvedarrow(theta,radius,center,varargin)
%CURVEDARROW plots a 2D curved arrow (surrogate to arrow.m)
% h=curvedarrow(theta [,radius,center,Length,BaseAngle,TipAngle,Width,Page,CrossDir])
%   theta = [thetastart thetastop] in degree
%   rradius = 1 (default)
%   center = [0 0] (default)
%   Length,BaseAngle,TipAngle,Width,Page,CrossDir: see arrowm
%   h = handles

%24/08/10 - INRA\Olivier Vitrac - rev. 

% argcheck
if nargin<1, theta = [180 100]; end
if nargin<2, radius = 1; end
if nargin<3, center = [0 0]; end
if nargin< 
if numel(theta)~=2, error('theta must be a 2x1 or 1x2 vector'), end
if numel(radius)~=1, error('radius must be a scalar'), end
if numel(center)~=2, error('center must be a 2x1 or 1x2 vector'), end
linewidth
m = 50;
theta = pi*theta/180;

% plots
THETA = linspace(theta(1),theta(2),m);
x = center(1)+cos(THETA);
y = center(2)+sin(THETA);
    