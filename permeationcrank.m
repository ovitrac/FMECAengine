function [flux_divl,Foout,flux_absout]=permeationcrank(Fo,n)
%PERMEATIONCRANK implements eq. 4.24a of Crank's book Mathematics of Diffusion, 2nd Ed., p 51 (+ eq. 4.23 p 50)
%   The main intent of this code is to offer high accurate solutions for validating numerical schemes
%   Assumptions: C1=1, C0=C2=0, Fo=Dt/l^2, Q = mass/m2 or mol/m2
%   Permeation Eq. Q/(lC1) = Fo -1/6 -2/pi^2*Sum(n=1:Inf,(-1)^n/n^2*exp(-n^2*pi^2*Fo)  
%   Sorption Eq. Mt/Moo = 1-8/pi^2 * Sum(n=0:Inf,1/(2*n+1)^2*exp(-(2*n+1)^2*pi^2*Fo) with Moo = l*((C1+C2)/2-C0)
%
% General syntax
%         flux_divl=permeationcrank()
%         flux_divl=permeationcrank(Fo [,n])
%         [flux_divl,Foout,flux_absout]=permeationcrank(Fo [,n])
%
%  Inputs (default values)
%         Fo: Fourier number (default = logspace(-6,1,1e4)')
%          n: number of terms in the serial expansion (default = 5000)
%
%  Outputs
% flux_divl: nFox1 vector of dimensionless permeated amount (eq 4.24a)
%        Fo: nFox1 vector Fourier numbers
%  flux_abs: nFox1 vector of dimensionless absorbed amount (eq. 4.23)
%
%       
% EXAMPLE to micmick the figure 4.4 of Crank (page 52)
%{
    [flux_divl,Fo,flux_abs]=permeationcrank;

    % permeate
    formatfig(figure,'paperposition',[ 4.8225   10.1726   11.3391    9.3323],'figname','permeationcranck_perm'), hold on
    plot(Fo,flux_divl,'b-','linewidth',2)
    plot(Fo,Fo-1/6,'r-','linewidth',1)
    formatax(gca,'xlim',[0 .6],'ylim',[0 .35],'fontsize',12)
    xlabel('Dt/l^2','fontsize',14)
    ylabel('Q_t/(lC_1)','fontsize',14)
    title('permeated amount','fontsize',14)

    % absorbate
    formatfig(figure,'paperposition',[ 4.8225   10.1726   11.3391    9.3323],'figname','permeationcranck_abs'), hold on
    plot(Fo,flux_abs,'b-','linewidth',2)
    formatax(gca,'xlim',[0 .6],'ylim',[0 1.07],'fontsize',12)
    xlabel('Dt/l^2','fontsize',14)
    ylabel('M_t/(l(C_1-C_2)/2-C_0)','fontsize',14)
    title('absorbed amount','fontsize',14)
    
    %print_pdf(600,get(gcf,'filename'),'','nocheck')
    %print_png(400,get(gcf,'filename'),'','',0,0,0)
%}

%INRA\Migration 2.0 - 01/12/12 - Olivier Vitrac - rev. 14/12:2017

% Revision history
% 14/12/2017 add absorption flux, updated example

% default values
n_default = 5e3;
Fo_default=[0;logspace(-6,1,1e4)'];

% arg check
if nargin<1, Fo=Fo_default; end
if nargin<2, n=n_default; end
Fo = Fo(:); nFo = length(Fo);

% do calc
[flux_divl,flux_abs] = deal(zeros(nFo,1));
for i=1:n
    flux_divl = flux_divl + ((-1)^i)/i^2 * exp(-Fo*i^2*pi^2);
    flux_abs = flux_abs + 1/(2*(i-1)+1)^2*exp(-Fo*(2*(i-1)+1)^2*pi^2);
end
flux_divl = Fo -1/6 -2/pi^2 * flux_divl;
flux_abs = 1 - 8/pi^2 * flux_abs;

% output
if nargout>1, Foout=Fo; end
if nargout>2, flux_absout = flux_abs; end
