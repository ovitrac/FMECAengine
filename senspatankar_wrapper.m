function S=senspatankar_wrapper(S0)
%SENSPATANKAR_WRAPPER to normalize data as senspatankar does
%   syntax: Sdimensionless=senspatankar_wrapper(Swithunits)
%
%   see also: SENSPATANKAR, SENSPATANKARC, SETOFFPATANKAR, SENSPATANKAR_RESTART


% MS-MATLAB-WEB 1.0 - 13/05/09 - Olivier Vitrac - rev. 24/10/2015

% revision history
% 28/04/11 add S.lengthscale
% 23/10/15 improved help
% 24/10/2015 add hasbeenwrapped (i.e. dimensionless)

S = S0;
[crit,iref] = min(S.D./(S.k.*S.l)); %#ok<ASGLU>
S.L = S.L * S.l(iref) / sum(S.l);
S.t = S.t*S.D(iref)/S.l(iref)^2;
S.timescale = S.l(iref)^2/S.D(iref);
S.lengthscale = S.l(iref); % reference length scale (added 28/04/11)
S.hasbeenwrapped = true;