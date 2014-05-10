function [D,alpha]=Dlimm(polymer,M,T)
% Dlimm returns the Limm's overestimate of diffusion coefficients (Limm et al., 1996)
%   syntax: D=Dlimm(polymer,M,T)
%           [D,alpha]=Dlimm(...) returns also the alpha values
%       polymer = 'LDPE'    'HDPE'    'PP'  
%       M = molecular mass
%       T = temperature in °C (default = 40°C)

% Migration 2.0 - 27/02/2013 - INRA\Olivier Vitrac - rev.

% Revision history

% definitions
data = struct(...
'LDPE'                   , struct('lnD0',4.16,'a',0.555,'b',1140.5),...
'HDPE'                   , struct('lnD0',0.9,'a',0.819,'b',1760.7),...
'PP'                     , struct('lnD0',-2.1,'a',0.597,'b',1335.7)...
    );

% arg check
if nargin<1, error('one argument is at least required'); end
if nargin<2, M = 100; end
if nargin<3, T = 40; end

% recursion
if iscell(polymer)
    D = [];
    for ip = 1:numel(polymer)
        Dtmp = Dpiringer(polymer{ip},M,T);
        D = [D;Dtmp(:)]; %#ok<AGROW>
    end
    return
end

% type control
if ischar(M), tmp = M; M=polymer; polymer = tmp; end
if ~ischar(polymer), error('polymer must a string or cell array of strings'), end

%polymer = upper(polymer); % removed  22/08/11
TK      = T+273.25;
lnD0 = data.(polymer).lnD0;
a = data.(polymer).a;
b = data.(polymer).b;
D       = 1e-4*exp(lnD0 + a*M.^(1/2) -b*M.^(1/3)./TK); % m2.s-1

% additional ouput
if nargout>1
    alpha =  sqrt(M).*a.*(-1.0./2.0)+(M.^(1.0./3.0).*b.*(1.0./3.0))./TK;
end
