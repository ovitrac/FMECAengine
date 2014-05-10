function [D,alpha]=Dhelmroth(polymer,M,T,proba)
% Dhelmroth returns the Limm's overestimate of diffusion coefficients (J Sci Food Agric 85:909–916 (2005))
%   syntax: D=Dhelmroth(polymer,M,T,proba)
%           alpha = Dhelmroth(...) returns also the alpha values
%       polymer = 'LDPE'    'HDPE'    'PP'  
%       M = molecular mass
%       T = temperature in °C (default = 40°C)

% Migration 2.0 - 27/02/2013 - INRA\Olivier Vitrac - rev.

% Revision history

% definitions
data = struct(...
'LDPE'                   , struct('a',1.2e-6,'b',0.37,'s',1.3),...
'HDPE'                   , struct('a',7.2e-7,'b',0.39,'s',1.6),...
'PP'                     , struct('a',1.9e-8,'b',0.36,'s',2)...
    );

% arg check
if nargin<1, error('one argument is at least required'); end
if nargin<2, M = 100; end
if nargin<3, T = 40; end
if nargin<4, proba = 0.5; end
if proba>1, proba = proba/100; end

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
if T~=23, D = M*NaN; return, end
a = data.(polymer).a;
b = data.(polymer).b;
s = data.(polymer).s;
M0 = 1;
D       = 1e-4*a*exp(-(M/M0).^b)*exp(s*norminv(1-(1-proba)/2,0,1)); % m2.s-1

% additional ouput
if nargout>1
    alpha =  (M.*b.*(M./M0).^(b-1.0))./M0;
end