function k = FMECAkair(temp)
%FMECAKAIR gives the theoretical kair value for FMECAengine() calculations
%   Syntax: k = FMECAkair(temp)
%           temp = temperature with units, such as '298 K', '25°C' (default = '40°C')
%   
%   Example:
%             k = FMECAkair('4°C')
%             k = FMECAkair('313 K')
%
%   See also: FMECAunit, FMECAvp, FMECADpiringer, FMECADfuller, FMECAgpolymer, FMECAKairP, load_chemspider, fmecaengine


% INRA\FMECAengine v 0.6 - 04/04/2015 - Olivier Vitrac - rev. 

% Revision history

% default
temp_default = '40°C';
R = 8.314472;

% arg check (basic one)
if nargin<1, temp = []; end
if isempty(temp), temp = temp_default; end

% set temperature in °C
T = FMECAunit('K',temp);

% k estimate
k = R*T;