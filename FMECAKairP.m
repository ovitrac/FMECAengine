function K = FMECAKairP(solute,temp,chi)
%FMECAAIRP calculates the theoretical air/polymer partition coefficient for FMECAengine() calculations
%   Syntax: K = FMECAKairP(solute,temp [,chi])
%           solute = string coding for a valid CAS number, chemical name, CSID, etc.
%           temp   = string coding for temperature including unit (default value = '25°C')
%           chi = Flory-Huggin coefficient in the polymer (default value = 0.27)
%           NOTE that solute and temp are STRINGS but that chi is a numeric value
%
%    Example:
%           K = FMECAKairP('anisole','-20°C');
%           K = FMECAKairP('decane','-20°C');
%           K = FMECAKairP('pentane'); % by default at 298 K
%           K = FMECAKairP('styrene','75°C');
%           K = FMECAKairP('BHT');
%
%   See also: FMECAunit, FMECADpiringer, FMECADfuller, FMECAkair, FMECAgpolymer, load_chemspider, fmecaengine


% INRA\FMECAengine v 0.6 - 04/04/2015 - Olivier Vitrac - rev. 14/04/2015


% Revision history
% 14/04/2015 fix proper output


% Default
temp_default = '25°C';

%arg check
if nargin<1, error('one argument is required'), end
if ~ischar(solute), error('mol must a be a string'), end
if nargin<2, temp = []; end
if isempty(temp), temp = temp_default; end
if ~ischar(temp), error('temp must be a string including the temperature value and its unit (e.g. 60°C)'), end
if nargin<3, chi = []; end

% check mol
dmol = load_chemspider(solute);
if isempty(dmol), error('the substance ''%s'' is not found',solute), end

% intermediate calculations
Vi = convertunit('cm3','m3',dmol.Properties.MolarVolume.value); % m3/mol
kair = FMECAkair(temp);       % RT
gpol =  FMECAgpolymer(chi);   % activity coefficient for the polymer
Pvsat = FMECAvp(solute,temp); % vapor pressure
K = Pvsat*Vi*gpol/kair;       % see PhD Thesis of Mai (2014), Eq. 2.25 (page 48)