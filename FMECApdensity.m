function d=FMECApdensity(polymer,temp)
%FMECApdensity returns the density in SI units (kg/m3) of commpn polymers
%   syntax: d=FMECApdensity(polymer [,temp])
%           polymer, temp = strings coding for the polymer name and temperature (with units)
%           by default, temp = '25°C'
%           FMECApdensity() returns the list of polymers
%
%   Examples:
%           d = FMECApdensity('PP')
%           d = FMECApdensity('PA 6,6','40°C')
%           d = FMECApdensity('PA 12','40°C')
%           d = FMECApdensity('ppvc','40°C')
%           d = FMECApdensity('PET>Tg','120°C')
%
%   Usage with fmecaengine (to calculate volume concentration in SI units)
%           FMECAunit('w','500 ppm')*FMECApdensity('HDPE','40°C')
%
%   List of polymers:
%         'LDPE'
%         'LLDPE'
%         'HDPE'
%         'PP'
%         'PPrubber'
%         'OPP'
%         'PS'
%         'HIPS'
%         'SBS'
%         'PET_sup_Tg'
%         'PET_inf_Tg'
%         'PBT'
%         'PEN'
%         'PA6'
%         'PA6_6'
%         'PA12'
%         'PVCrigid'
%         'pPVC'
%         'AdhesiveNaturalRubber'
%         'AdhesiveSyntheticRubber'
%         'AdhesiveEVA'
%         'AdhesiveVAE'
%         'AdhesivePVAC'
%         'AdhesiveAcrylate'
%         'AdhesivePU'
%         'Paper'
%         'Cardboard_polarmigrant'
%         'Cardboard_apolarmigrant'
% 
%
%   See also: FMECAunit, FMECADpiringer, FMECADfuller, FMECAkair, FMECAgpolymer, FMECAKairP, load_chemspider, fmecaengine
%
%   Source of data: Referex\Chemistry\Polymer data handbook 1999 - Mark.pdf

% INRA\FMECAengine v 0.6 - 04/04/2015 - Olivier Vitrac - rev. 

% Revision history

% default
temp_default = '25°C';
replacements = {...
    '>Tg' '_sup_Tg'
    'overTg' '_sup_Tg'
    '<Tg' '_below_Tg'
    '[\s\_\-\,\;]',''
};

% definitions
data = struct(...
'LDPE'                   , struct('dmin',0.910,'dmax',0.935,  'z',30e-5),...
'LLDPE'                  , struct('dmin',0.863,'dmax',0.885,  'z',60e-5),...
'HDPE'                   , struct('dmin',0.94,'dmax',0.96,    'z',10e-5),... @(tC) 0.0894/(4578-22.7*tC)
'PP'                     , struct('dmin',0.90,'dmax',0.91,    'z',66e-5),...
'PPrubber'               , struct('dmin',0.85,'dmax',0.86,    'z',90e-5),...
'OPP'                    , struct('dmin',0.93,'dmax',0.946,   'z',20e-5),...
'PS'                     , struct('dmin',1.04,'dmax',1.065,   'z',10e-5),...
'HIPS'                   , struct('dmin',1.111,'dmax',1.1127, 'z',4e-5),...
'SBS'                    , struct('dmin',0.92,'dmax',0.95,    'z',10e-5),...
'PET_sup_Tg'             , struct('dmin',1.40 ,'dmax',1.40,   'z',60e-5),...
'PET_inf_Tg'             , struct('dmin',1.41 ,'dmax',1.41,   'z',10e-5),...
'PBT'                    , struct('dmin',1.33 ,'dmax',1.34,   'z',10e-5),...
'PEN'                    , struct('dmin',1.34 ,'dmax',1.3471, 'z',10e-5),...
'PA6'                    , struct('dmin',1.13 ,'dmax',1.13,   'z', 83e-5),...
'PA6_6'                  , struct('dmin',1.27 ,'dmax',1.27,   'z', 30e-5),...
'PA12'                   , struct('dmin',1.01 ,'dmax',1.03,   'z', 20e-5),...
'PVCrigid'               , struct('dmin',1.36, 'dmax',1.35,   'z', 45e-5),...
'pPVC'                   , struct('dmin',1.10, 'dmax',1.25,   'z', 90e-5),...
'AdhesiveNaturalRubber'  , struct('dmin',0.91, 'dmax',0.93,   'z',60e-5),...
'AdhesiveSyntheticRubber', struct('dmin',0.92, 'dmax',0.94,    'z',60e-5),...
'AdhesiveEVA'            , struct('dmin',0.93 ,'dmax',0.94,    'z', 40e-5),...
'AdhesiveVAE'            , struct('dmin',0.94 ,'dmax',0.94,    'z', 40e-5),...
'AdhesivePVAC'           , struct('dmin',0.94 ,'dmax',0.94,    'z', 40e-5),...
'AdhesiveAcrylate'       , struct('dmin',1.22 ,'dmax',1.22,    'z', 12e-5),...
'AdhesivePU'             , struct('dmin',0.98 ,'dmax',1.24,    'z',60e-5),...
'Paper'                  , struct('dmin',0.3  ,'dmax',0.4,'z',0),...
'Cardboard_polarmigrant' , struct('dmin',0.4  ,'dmax',0.4,'z',0),... polar substances
'Cardboard_apolarmigrant', struct('dmin',0.4 ,'dmax',0.3,'z',0) ... hydrophobic substances
    );

% arg check
if nargin<1 && ~nargout, d = fieldnames(data); return, end
if nargin<1, error('one argument is at least required'); end
if nargin<2, temp = ''; end
if isempty(temp), temp = temp_default; end
if ~ischar(polymer), error('polymer must be a string'), end
if ~ischar(temp), error('temperature must be a string including units'), end

% polymer hash and polymer identification
polymerlist = fieldnames(data);
polymerhash = lower(regexprep(polymerlist,'_',''));
polpattern = lower(regexprep(polymer,replacements(:,1),replacements(:,2)));
ifound = find(ismember(polymerhash,polpattern));
if length(ifound)~=1
    dispf('FMECApdensity:: use approximative search for ''%s''',polymer)
    ifound = find(~cellfun(@isempty,uncell(regexp(polymerhash,['^' polpattern]))));
    foundpol = polymerlist(ifound);
    nfound = length(ifound);
    if nfound==0, error('no polymer/material found for ''%s''',polymer), end
    if nfound>1
        dispf('\nFMECApdensity(''%s'',''%s'')\n\t%d polymer/materials found:',polymer,temp,nfound)
        dispf('\t>\t%s\n',foundpol{:})
        error('%d polymer/materials matches found for the estimation of density of ''%s''',nfound,polymer)
    end
end

% temperature
T0 = 25;
T = FMECAunit('C',temp);
dref = data.(polymerlist{ifound});
d = 1/(mean([1/dref.dmin,1/dref.dmax])*(1+(T-T0)*dref.z));

% final conversion
d = convertunit('g/cm3','kg/m3',d);
