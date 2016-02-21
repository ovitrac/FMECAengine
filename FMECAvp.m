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
%           Pvsat = FMECAvp('water','100°C')
%
%   See also: FMECAunit, FMECADpiringer, FMECADfuller, FMECAKairP, load_chemspider, fmecaengine


% INRA\FMECAengine v 0.6 - 04/04/2015 - Olivier Vitrac - rev.  08/02/2016

% Revision history
% 07/02/2016 add UserProperties for VP and MP if missing (now FMECAKairP limonene works)
% 07/02/2016 intercept errors in ModifiedGrainMethod(), FMECAKairP ethane displays a warning
% 08/02/2016 check the version of Chemspider


% Default
Tdefault = '25°C';
tempclean = @(Tstr) regexprep(Tstr,{'\s' '°C'},{'' 'degC'});

%arg check
if nargin<1, error('one argument is required'), end
if ~ischar(mol), error('mol must a be a string'), end
if nargin<2, T = []; end
if isempty(T), T = Tdefault; end
if ~ischar(T), error('T must be a string including the temperature value and its unit (e.g. 60°C)'), end

% check mol
dmol = load_chemspider(mol);
if isempty(dmol), error('the substance ''%s'' is not found',mol), end
if ~isfield(dmol,'EPI'), error('EPI data are missing, check your version of load_chemspider'), end
if isempty(dmol.EPI) || ~isfield(dmol.EPI,'BoilingPt_MeltingPt_VaporPressureEstimations'), error('no EPI data for the substance ''%s''',dmol.Name), end
thermo = dmol.EPI.BoilingPt_MeltingPt_VaporPressureEstimations;

% reference pressure in units required by ModifiedGrainMethod()
VP = thermo.VP; % reference pressure
VP.interpreted = regexp(VP.unit,'(?<Punit>[\w\s]*)\s*at\s*(?<T>\d+)\s*(?<Tunit>.*)','names');
P0 = convertunit(regexprep(VP.interpreted.Punit,'\s',''),... initial unit
                 'atm',... destination unit
                 VP.value); % reference P value
T0 = convertunit(tempclean(VP.interpreted.Tunit),...
                 'K',...
                 str2double(VP.interpreted.T),...
                 'abstemp');
             
% reference boiling point in units required by ModifiedGrainMethod()
if isfield(thermo,'BP') && ~isempty(thermo.BP)
    BP = thermo.BP;
elseif isfield(dmol.UserProperties,'ExperimentalBoilingPoint') && ~isempty(dmol.UserProperties.ExperimentalBoilingPoint)
    BP = struct('value',mean(dmol.UserProperties.ExperimentalBoilingPoint.value),...
                'unit',dmol.UserProperties.ExperimentalBoilingPoint.unit{1});
else
    error('FMECAvp:: no boiling point (BP) for ''%s'' in EPIsuite and in UserProperties',mol)
end
Tb = convertunit(tempclean(BP.unit),...
                 'K',...
                 BP.value,...
                 'abstemp');

% reference melting point in units required by ModifiedGrainMethod()
if isfield(thermo,'MP') && ~isempty(thermo.MP)
    MP = thermo.MP;
elseif isfield(dmol.UserProperties,'ExperimentalMeltingPoint')
    dispf('FMECAvp:: no melting point (MP) for ''%s'' in EPIsuite, switch to user properties',mol)
    MP = struct('value',mean(dmol.UserProperties.ExperimentalMeltingPoint.value),...
                'unit',dmol.UserProperties.ExperimentalMeltingPoint.unit{1});
else
    error('FMECAvp:: no melting point for ''%s'' in EPIsuite and in UserProperties',mol)
end
Tm = convertunit(tempclean(MP.unit),...
                 'K',...
                 MP.value,...
                 'abstemp');

% pressure calculations at temperature T
TK = FMECAunit('K',T);
Pvsat = convertunit('atm','Pa',ModifiedGrainMethod(P0,T0,TK,Tb,Tm));

% detection of complex pressure
if ~isreal(Pvsat)
    dispf('FMECAcp:: ModifiedGrainMethod() fails, a complex Pvsat has been detected for ''%s''',mol)
    Pvsat = real(Pvsat);
end
