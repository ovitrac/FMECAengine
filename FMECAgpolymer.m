function k = FMECAgpolymer(chi)
%FMECAKAIR gives the theoretical gpolymer (activity coefficient) value for FMECAengine() calculations
%   Syntax: k = FMECAgpolymer(chi)
%           chi = numerical chi value (default = 0.27)
%   
%   Example:
%             g = FMECAgpolymer
%             g = FMECAgpolymer(2)
%
%
%   See also: FMECAunit, FMECAvp, FMECADpiringer, FMECADfuller, FMECAkair, FMECAKairP, load_chemspider, fmecaengine


% INRA\FMECAengine v 0.6 - 04/04/2015 - Olivier Vitrac - rev. 

% Revision history

% default
chi_default = 0.27;

% arg check (basic one)
if nargin<1, chi = []; end
if isempty(chi), chi = chi_default; end

% G estimate
k = exp(1+chi);