function ye = sumexpval(R,t)
% SUMEXPVAL évalue une somme d'exponentielle décroissante (avec biais éventuel)
%		ye = sumexpval(R,t)

% Woodox 1.21 - 03/04/01 - Olivier Vitrac - rev. 28/04/03

% Revision history
% 28/04/04 arg check added, add type 'incexp'

% arg check
if nargchk(2,2,nargin), error('2 arguments are required'), end
if ~isstruct(R) | ~isfield(R,'type') ~isfield(R,'Cend'), error('R must be a structure created by sumexpfit'), end
t = t(:);

% calculations
m = length(t);
n = R.order; %length(R.l);
LAMBDAS = zeros(m,n);
for i = 1:n, LAMBDAS(:,i) = exp(-R.l(i)*t); end
ye = LAMBDAS*R.C(1:n)' + R.Cinf;
if strcmp(R.type,'incexp')
    ye = R.Cend - ye;
end
