function S=senspatankar_wrapper(S0)
%SENSPATANKAR_WRAPPER to normalize data as senspatankar does

% MS-MATLAB-WEB 1.0 - 13/05/09 - Olivier Vitrac - rev. 28/04/11

% revision history
% 28/04/11 add S.lengthscale

S = S0;
[crit,iref] = min(S.D./(S.k.*S.l)); %#ok<ASGLU>
S.L = S.L * S.l(iref) / sum(S.l);
S.t = S.t*S.D(iref)/S.l(iref)^2;
S.lengthscale = S.l(iref); % reference length scale (added 28/04/11)