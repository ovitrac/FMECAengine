function Pvsat = FMECAvp(mol,T)
%FMECAVP vapor pressure calculation based on ModifiedGrainMethod(P0,T0,T,Tb,Tm) and EPISUITE/MPBPWIN
%   Syntax: Pvsat = FMECAvp(mol,T)
%           mol = string coding for a valid CAS number, chemical name, CSID, etc.
%           T   = string coding for temperature including unit (default value = '25°C')
%           Pvsat = liquid vapor pressure (numeric value in Pa)
%
%    Example:
%           Pvsat = FMECAvp('anisole','-20°C');
%           Pvsat = FMECAvp('decane','-20°C');
%           Pvsat = FMECAvp('pentane'); % by default at 298 K
%           Pvsat = FMECAvp('styrene','75°C');
%           Pvsat = FMECAvp('BHT');
%
%   See also: FMECAunit, FMECADpiringer, FMECADfuller, FMECAKairP, load_chemspider, fmecaengine


% INRA\FMECAengine v 0.6 - 04/04/2015 - Olivier Vitrac - rev. 

% Revision history


% Default
Tdefault = '25°C';

%arg check
if nargin<1, error('one argument is required'), end
if ~ischar(mol), error('mol must a be a string'), end
if nargin<2, T = []; end
if isempty(T), T = Tdefault; end
if ~ischar(T), error('T must be a string including the temperature value and its unit (e.g. 60°C)'), end

% check mol
dmol = load_chemspider(mol);
if isempty(dmol), error('the substance ''%s'' is not found',mol), end
if isempty(dmol.EPI), error('no EPI data for the substance ''%s''',dmol.Name), end
thermo = dmol.EPI.BoilingPt_MeltingPt_VaporPressureEstimations;

% reference pressure in units required by ModifiedGrainMethod()
VP = thermo.VP; % reference pressure
VP.interpreted = regexp(VP.unit,'(?<Punit>[\w\s]*)\s*at\s*(?<T>\d+)\s*(?<Tunit>.*)','names');
P0 = convertunit(regexprep(VP.interpreted.Punit,'\s',''),... initial unit
                 'atm',... destination unit
                 VP.value); % reference P value
T0 = convertunit(regexprep(VP.interpreted.Tunit,'\s',''),...
                 'K',...
                 str2double(VP.interpreted.T),...
                 'abstemp');
             
% reference boiling point in units required by ModifiedGrainMethod()
BP = thermo.BP;
Tb = convertunit(regexprep(BP.unit,'\s',''),...
                 'K',...
                 BP.value,...
                 'abstemp');

% reference melting point in units required by ModifiedGrainMethod()
MP = thermo.MP;
Tm = convertunit(regexprep(MP.unit,'\s',''),...
                 'K',...
                 MP.value,...
                 'abstemp');

% pressure calculations at temperature T
TK = FMECAunit('K',T);
Pvsat = convertunit('atm','Pa',ModifiedGrainMethod(P0,T0,TK,Tb,Tm));
