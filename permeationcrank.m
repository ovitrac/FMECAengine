function [flux_divl,Foout]=permeationcrank(Fo,n)
%PERMEATIONCRANK implements eq. 4.24a of Crank's book Mathematics of Diffusion, 2nd Ed., p 51
%   The main intent of this code is to offer high accurate solutions for validating numerical schemes
%   Assumptions: C1=1, C0=C2=0, Fo=Dt/l^2, Q = mass/m2 or mol/m2
%   Original Eq. Q/(lC1) = Fo -1/6 -2/pi^2*Sum(n=1:Inf,(-1)^n/n^2*exp(-n^2*pi^2*Fo)  
% syntax: flux_divl=permeationcrank;
%         flux_divl=permeationcrank(Fo);
%         flux_divl=permeationcrank(Fo,n);
%         flux_divl nFox1 vector of dimensionless flux density x l
%   default values
%       Fo=logspace(-6,1,1e4)';
%       n=5000
%       
% EXAMPLE to micmick the figure 4.4 of Crank (page 52)
% [flux_divl,Fo]=permeationcrank;
% figure, hold on
% plot(Fo,flux_divl,'b-','linewidth',2)
% plot(Fo,Fo-1/6,'r-','linewidth',1)
% formatax(gca,'xlim',[0 .6],'ylim',[0 .35],'fontsize',12)
% xlabel('Dt/l^2','fontsize',14)
% ylabel('Q_t/(lC_1)','fontsize',14)
% definitions
% set(gcf,'paperposition',[ 4.8225   10.1726   11.3391    9.3323])
% print_pdf(600,'permeation','','nocheck')
% print_png(400,'permeation','','',0,0,0)

%INRA\Migration 2.0 - 01/12/12 - Olivier Vitrac - rev.

% default values
n_default = 5e3;
Fo_default=logspace(-6,1,1e4)';

% arg check
if nargin<1, Fo=Fo_default; end
if nargin<2, n=n_default; end
Fo = Fo(:); nFo = length(Fo);

% do calc
flux_divl = zeros(nFo,1);
for i=1:n
    flux_divl = flux_divl + ((-1)^i)/i^2 * exp(-Fo*i^2*pi^2);
end
flux_divl = Fo -1/6 -2/pi^2 * flux_divl;

% output
if nargout>1, Foout=Fo; end
