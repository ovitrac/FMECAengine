function Dair = Dfuller(formula,T,P,A,Hc)
% DFULLER returns the diffusion coefficients of gases in air (Fuller et al. (1965, 1966,1969))
% DAB = (0.00143*T^1.75)/(P*sqrt(MAB)*(vA^(1/3)+vB^(1/3))^2
% MAB = 2/(1/MA+1/MB)
%
%   syntax: D = Dfuller(gas,formula,T,P) in m/s2
% formula = formule of soulte
%       T = temperature in °C (default = 25)
%       P = pressure in Pa (default 1e5 Pa)
%
%   Source: The Properties of Gases and Liquids, Fifth Edition Bruce E. Poling, John M. Prausnitz, John P. O’Connell
%
% Migration 2.0 - 31/03/2014 - INRA\Olivier Vitrac - rev. 

%% argcheck
if nargin < 1, error('1 argument is required'), end
if nargin < 2 || isempty(T), T = 25; end
if nargin < 3 || isempty(P), P = 1e5 ; end
if nargin < 4 || isempty(A), A = 6*12+5; end
if nargin < 5 || isempty(Hc), Hc = NaN; end


%% definitions
Vair = 19.7;
Mair = 0.8*14*2+0.2*16*2;
AtomicDiffusionVolumes = {...
    'C'    'H'   'O'   'N'    'A'   'Hc'    'F'   'Cl'  'Br'  'I'   'S'
    15.9   2.31  6.11  4.54  -18.3  -18.3  14.7  21    21.9  29.8  22.9 }';
MolecularMass = {...
    'C'       'H'      'O'       'N'       'F'         'Cl'         'Br'  'I'   'S'  'A' 'Hc'
     12.011  1.008     15.999    14.007    19.997999    35.452999    80   126   32.06  A  Hc}';
[~,isort] = sort(cellfun(@length,AtomicDiffusionVolumes(:,1)),'descend');
AtomicDiffusionVolumes = AtomicDiffusionVolumes(isort,:);
patterns = cellfun(@(s) sprintf('(?<!#)(%s\\d*)(?!#)',s),AtomicDiffusionVolumes(:,1),'UniformOutput',false);
Vcode = cell2struct(AtomicDiffusionVolumes(:,2),AtomicDiffusionVolumes(:,1),1);
Mcode = cell2struct(MolecularMass(:,2),MolecularMass(:,1),1);

% interpreter
protected = regexprep(regexprep(formula,patterns,' #$1# '),'([A-Za-z])#','$1\1#');
code = regexp(protected,'\#(?<atom>[A-Za-z]+)(?<number>\d+)','names');
N = cellfun(@str2double,{code.number});
Va = cellfun(@(c) Vcode.(c),{code.atom});
Ma = cellfun(@(c) Mcode.(c),{code.atom});

% do calculation
Vsolute = N*Va';
Msolute = N*Ma';
MAB = 2/(1/Mair+1/Msolute);
Dair = 1e-4*(0.00143.*(T+273.15).^1.75)./(1e-5*P*sqrt(MAB)*(Vair^(1/3)+Vsolute^(1/3))^2);