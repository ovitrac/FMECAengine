function Sr = senspatankar_restart(S,R0)
%SENSPATANKAR_RESTART prepares simulation inputs to reuse previous simulation results
%  SYNTAX: srestart = senspatankar_restart(s,r0)
%   USAGE:  r  = senspatankarC(senspatankar_wrapper(senspatankar_restart(s,r0)))
%
%   see also: SENSPATANKAR, SENSPATANKARC, SETOFFPATANKAR, SENSPATANKAR_WRAPPER
%
% =================================================================================
% | EXAMPLE 1 |       2 layers, FB (contact) layer removed
% =================================================================================
%{
s0 = struct('Bi',1e3,...
            'k',[1 1],...
            'D',[1e-16 1e-13],...
            'k0',1,...
            'l', [10e-6 200e-6],...
            'L', 0.1,...
            'C0',[0  1000],...
            't', linspace(0,100*24*3600,100));
s = s0; s.C0 = [-1 0];                    % remove layer 1
s.t=linspace(100*24*3600,120*24*3600,20) % continuous time scale
r0 = senspatankar(senspatankar_wrapper(s0));
r  = senspatankarC(senspatankar_wrapper(senspatankar_restart(s,r0)));
figure % -------------- plot concentration profiles
subplot(211), plot(r0.x*r0.F.lengthscale,r0.Cx(1:10:end,:)','b-'), title('EX1/step1: every 10 days'), ylabel('profile')
subplot(212), plot(r.x*r.F.lengthscale,r.Cx','r-'), title('EX1/step2: every day'), xlabel('position'), ylabel('profile')
figure % -------------- plot kinetics
plot(r0.t*r0.timebase/(24*3600),r0.CF,'b-',r.t*r.timebase/(24*3600),r.CF,'r-'), title('Example 1'), xlabel('time'), ylabel('C_F')
%}
% =================================================================================
% | EXAMPLE 2 |  3 layers (symmetric assembling ABC, source B is removed at step 2)
% =================================================================================
%{
s0 = struct('Bi',1e3,...
            'k',[1 1 1],...
            'D',[1e-16 1e-13 1e-16],...
            'k0',1,...
            'l', [10e-6 200e-6 10e-6],...
            'L', 0.1,...
            'C0',[0  1000 0],...
            't', linspace(0,100*24*3600,100));
s = s0; s.C0 = [0 -1 0];                    % remove source (layer 2)
s.t=linspace(100*24*3600,150*24*3600,50) % continuous time scale
r0 = senspatankar(senspatankar_wrapper(s0));
r  = senspatankarC(senspatankar_wrapper(senspatankar_restart(s,r0)));
figure % -------------- plot concentration profiles
subplot(211), plot(r0.x*r0.F.lengthscale,r0.Cx(1:2:end,:)','b-'), title('EX2/step1: every 2 day'), ylabel('profile')
subplot(212), plot(r.x*r.F.lengthscale,r.Cx','r-'), title('EX2/step2: every day'), xlabel('position'), ylabel('profile')
figure % -------------- plot kinetics
plot(r0.t*r0.timebase/(24*3600),r0.CF,'b-',r.t*r.timebase/(24*3600),r.CF,'r-'), title('Example 2'), xlabel('time'), ylabel('C_F')
%}
% =================================================================================
% | EXAMPLE 3 |        3 layers with reservoir layer removed at step 2
% =================================================================================
%{
s0 = struct('Bi',1e3,...
            'k',[1 1 1],...
            'D',[1e-16 1e-16 1e-13],...
            'k0',1,...
            'l', [10e-6 10e-6 200e-6],...
            'L', 0.1,...
            'C0',[0 0 1000],...
            't', linspace(0,100*24*3600,100));
s = s0; s.C0 = [0 0 -1];                    % remove source (layer 2)
s.t=linspace(100*24*3600,300*24*3600,200) % continuous time scale
r0 = senspatankar(senspatankar_wrapper(s0));
r  = senspatankarC(senspatankar_wrapper(senspatankar_restart(s,r0)));
figure % -------------- plot concentration profiles
subplot(211), plot(r0.x*r0.F.lengthscale,r0.Cx(1:2:end,:)','b-'), title('EX3/step1: every 2 day'), ylabel('profile')
subplot(212), plot(r.x*r.F.lengthscale,r.Cx(1:10:end,:)','r-'), title('EX3/step2: every 10 days'), xlabel('position'), ylabel('profile')
figure % -------------- plot kinetics
plot(r0.t*r0.timebase/(24*3600),r0.CF,'b-',r.t*r.timebase/(24*3600),r.CF,'r-'), title('Example 3'), xlabel('time'), ylabel('C_F')
%}
%

% MS-MATLAB-WEB 2.0 - 23/10/15 - Olivier Vitrac - rev. 24/10/2015

% revision history
% 23/10/2015 release candidate 
% 24/10/2015 check hasbeenwrapped, propagate lengthscale
% 24/10/2015 add example 3 to verify the continuity of fluxes with time

% Default
interpmethod = 'pchip';

% check that S has not been wrapped
if isfield(S,'hasbeenwrapped') && S.hasbeenwrapped
    error('SENSPATANKAR_WRAPPER has been called before SENSPATANKAR_RESTART (this function), please do the reverse.')
end
if ~isfield(R0.F,'hasbeenwrapped') || ~R0.F.hasbeenwrapped
    error('R0 should have been called through SENSPATANKAR_WRAPPER')
end
    
% Extension
Sr  = S;
Sr.restart = struct(...
        'x',R0.x,...
        'C',interp1(sqrt(R0.tC*R0.timebase),R0.Cx,sqrt(S.t(1)),interpmethod),...
        'CF',interp1(R0.t*R0.timebase,R0.CF,S.t(1),interpmethod),...
        'layerid',R0.layerid,...
        'xlayerid',R0.xlayerid,...
        'C0eq',R0.C0eq,...
        'lengthscale',R0.F.lengthscale ...
                );