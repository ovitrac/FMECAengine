function D = FMECADpiringer(polymer,solute,temp)
%FMECADPIRINGER wrapper to Dpringer() for FMECAengine() calculations
%   Syntax: D = FMECADpiringer(polymer,solute [,temp])
%           polymer, solute, temp are all STRINGS
%           polymer = polymer code (see Diringer for calculation)
%           solute = any valid CAS, CSID, chemical name, InChIKey, etc..
%           temp = temperature with units, such as '298 K', '25°C' (default = '40°C')
%   
%   Example:
%             D = FMECADpiringer('PP','anisole','4°C')
%             D = FMECADpiringer('PP','decane','313 K')
%             D = FMECADpiringer('LDPE','BHT','313 K')
%             D = FMECADpiringer('PA 6,6','anisole','4°C')
%
%   See also: FMECAunit, FMECAvp, FMECADfuller, FMECAKairP, load_chemspider, fmecaengine


% INRA\FMECAengine v 0.6 - 04/04/2015 - Olivier Vitrac - rev. 05/04/2015

% Revision history
% 05/04/2015 check polymer name based on the definitions in Dpiringer (original code from FMECApdensity)

% default
temp_default = '40°C';
replacements = {...
    '>Tg' '_sup_Tg'
    'overTg' '_sup_Tg'
    '<Tg' '_below_Tg'
    '[\s\_\-\,\;]',''
};

% arg check (basic one)
if nargin<2, error('two arguments are at least required'), end
if ~ischar(polymer), error('polymer must be a string'), end
if ~ischar(solute), error('solute must be a string'), end
if nargin<3, temp = []; end
if isempty(temp), temp = temp_default; end

% check polymer name
polymerlist = Dpiringer();
polymerhash = lower(regexprep(polymerlist,'_',''));
polpattern = lower(regexprep(polymer,replacements(:,1),replacements(:,2)));
ifound = find(ismember(polymerhash,polpattern));
if length(ifound)~=1
    dispf('FMECADpiringer:: use approximative search for ''%s''',polymer)
    ifound = find(~cellfun(@isempty,uncell(regexp(polymerhash,['^' polpattern]))));
    foundpol = polymerlist(ifound);
    nfound = length(ifound);
    if nfound==0, error('no polymer/material found for ''%s''',polymer), end
    if nfound>1
        dispf('\nFMECADpiringer(''%s'',''%s'',''%s'')\n\t%d polymer/materials found:',polymer,solute,temp,nfound)
        dispf('\t>\t%s\n',foundpol{:})
        error('%d polymer/materials matches found for the overestimation of D in ''%s''',nfound,polymer)
    end
end
polymer = polymerlist{ifound};

% retrieve the substance and M
dmol = load_chemspider(solute);
if isempty(dmol), error('the substance ''%s'' is not found',solute), end
M = str2double(dmol.QuickMass.molecularmass);

% set temperature in °C
T = FMECAunit('C',temp);

% D "overestimate"
D = Dpiringer(polymer,M,T);
