function D = FMECADfuller(solute,temp)
%FMECADFULLER wrapper to Dfuller() for FMECAengine() calculations
%   Syntax: D = FMECADfuller(polymer,solute [,temp])
%           solute, temp are all STRINGS
%           solute = any valid CAS, CSID, chemical name, InChIKey, etc..
%           temp = temperature with units, such as '298 K', '25°C' (default = '40°C')
%   
%   Example:
%             D = FMECADfuller('anisole','4°C')
%             D = FMECADfuller('decane','313 K')
%             D = FMECADfuller('BHT','313 K')
%
%   See also: FMECAunit, FMECAvp, FMECADpiringer, FMECAKairP, load_chemspider, fmecaengine


% INRA\FMECAengine v 0.6 - 04/04/2015 - Olivier Vitrac - rev. 

% Revision history

% default
temp_default = '40°C';
P = 1e5; % Pa

% arg check (basic one)
if nargin<1, error('one argument is at least required'), end
if ~ischar(solute), error('solute must be a string'), end
if nargin<2, temp = []; end
if isempty(temp), temp = temp_default; end

% retrieve the substance and M
dmol = load_chemspider(solute);
if isempty(dmol), error('the substance ''%s'' is not found',solute), end
formula = dmol.QuickMass.formula;

% set temperature in °C
T = FMECAunit('C',temp);

% D estimate at 1e5 Pa
D = Dfuller(formula,T,P);