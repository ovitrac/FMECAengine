%i;exec('matlab2scilab.m',-1)
function [resode,respde,resanalytic]=sensc(methode,F,ploton,dispon)
%  // SENSC simule un transfert RD+RK+RH lors d'un contact solide fluide avec des isothermes arbitraires (BET, GP...)
%  //   Syntaxe: F = sensc()
%  //   Syntaxe: [RESODE,RESPDE,RESANALYTIC] = SENSC(MODE,F[,PLOTON,DISPON])
%
% Values for MODE :
%     'init'     : constructor of F (default value of MODE).
%     'auto'     : tout.
%     'ode'      : Algebraic solution based on Vitrac Hayert 2006. Result in resode. Attention, ne marche pas !
%     'pde'      : Finite element solution. Result in respde.
%     'analytic' : Analytic solution based on Vitrac Goujot 2008 (in preparation). Result in resanalytic. If first/second letter is uppercase, use only first/second matrix product at corners. If third letter is uppercase of if ploton is true, combine F.t and F.x in resanalytic.tC
%     'odepde'   : ode and pde.
%     'all'      : ode and pde and analytic.
%
% Values for F :
%     sol : le solide charge initialement avec la concentration C0sol mol/m^3 initiale.
%     flu : le fluide charge initialement avec la concentration C0flu mol/m^3 initiale. TODO BUG: Tres mal gere si non nul.
%
% Values for ploton :
%     false      : no figure (default).
%     true       : various plots for visual check.
%     2          : various plots for visual check, do not close or clear figures.
%
% Values for dispon :
%     false      : no comment (default).
%     true       : various output showing progress.
%
% Debug is triggered if some letters of 'analytic' are upper case, with following meaning:
%     'Analytic' : include alternate matrix product.
%     'ANalytic' : exclude stable matrix product.
%     'anAlytic' : include profiles for all user timesteps F.t and some additional times around corners (automatic if ploton is set).
%     'anaLytic' : use Sagiv product with inverse matrix.
%     'anaLYtic' : compare Sagiv with stable matrix product.
%     'anaLyTic' : Use Sagiv enhanced with linear inversion.
%     'analytIc' : resanalytic.err gives clues about precision near t=t0, but no more around corners in profiles.
%     'anaLytiC' : use Sagiv product, with inv instead of \ ; also give plot debug about eigenvalue computatoin;
% SENS 1.0 - 22/02/07 - Olivier Vitrac - rev.
% SENSC 1.1 - 09/05/07 - Daniel Goujot - rev.
% SENSC 1.2 - 11/10/07 - Daniel Goujot - rev.
%
% revision history
% SENSC 1.1 - 09/05/07 - isotherme a coin.
% SENSC 1.2 - 11/10/07 - solution analytique.
%i;// a utiliser avec .PvsatsurRT
%i;// appel sous matlab : cd K:/recherCHE/matlab/senscrep,F=sensc,betiso=feval(F.iso.helperfunctions.baldev02isotherms,'BET',feval(F.iso.helperfunctions.baldev02STARCHandPE,'BET',.4));betiso.awminawmax(3)=-2-1/10;F.iso.S=feval(F.iso.helperfunctions.approximate_isotherm,betiso);F.C0sol=5;F.L=.1;[resode,respde,resanalytic]=sensc('anaLYTic',F,true,true) ; OK, marche 20080611
%i;// 'Pvsat(273.15+27) et 1atm en pascal :';pvsat27=3564.9138;atm=101338.52;
%i;// density of LDPE http://www.freepatentsonline.com/EP0152865A1.html 922 kg/m^3.
%i;// wikipedia density air:  density of air at 1atm at 20°C under 1atm:  1.2041 kg/m^3
%i;// wikipedia density air:  \rho~_{_{humid~air}} = \frac{(p-aw*pvsat)}{RsurMair \cdot T} + \frac{aw*pvsat}{RsurMeau \cdot T}  en kg/m^3. avec RsurMair=287.05 et RsurMeau=461.495.
%i;// /usr/bin/bc : 101338.52/287.05/(273.15+20) donne 1.20427898432070391574
%i;// /usr/bin/bc : 101338.52/287.05/(273.15+27) donne 1.17619318425325454905
%i;// /usr/bin/bc : (101338.52-.2*3564.9)/287.05/(273.15+27)+.2*3564.9/461.495/(273.15+27) donne 1.17306513937994601324
%i;// 
%i;// http://www.gpa.uq.edu.au/UQweather/airdensity.htm :  Air Density (in Kg/m3) = 1.2929 X 273.13 X ( AP - ( SVP x RH )) ( T + 273.13) 760
%i;// 
%i;// appel sous matlab : cd K:/recherCHE/matlab/senscrep,F=sensc,betiso=feval(F.iso.helperfunctions.baldev02isotherms,'BET',feval(F.iso.helperfunctions.baldev02STARCHandPE,'BET',.2));betiso.awminawmax(1)=0;betiso.awminawmax(3)=-1/10;F.iso.S=feval(F.iso.helperfunctions.approximate_isotherm,betiso);F.C0sol=1;F.L=.01;F.Bi=10;F.t=[1 10]';[resode,respde,resanalytic]=sensc('anaLYTic',F,true,true)
%i;// figure(1);clf;hold on;F=sensc;for i=.1:.1:.5;betiso=feval(F.iso.helperfunctions.baldev02isotherms,'BET',feval(F.iso.helperfunctions.baldev02STARCHandPE,'BET',i));betiso.awminawmax(3)=-1/10;F.iso.S=feval(F.iso.helperfunctions.approximate_isotherm,betiso);plot([-.1,F.iso.S.aws],[i,F.iso.S.Cs/F.iso.S.Cs(end),i]);end;
%i;// figure(2);clf;hold on;densiteLDPE=922;pvsat27=3564.9138;RsurMeau=461.495;F=sensc;for i=.1:.1:.5;betiso=feval(F.iso.helperfunctions.baldev02isotherms,'BET',feval(F.iso.helperfunctions.baldev02STARCHandPE,'BET',i));betiso.awminawmax(3)=-1/10;F.iso.S=feval(F.iso.helperfunctions.approximate_isotherm,betiso);plot(F.iso.S.aws,[pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K/densiteLDPE,i/100]);end;
%i;// clear;figure(10);clf;liste={'bet' 'smith' 'halsey' 'caurie' 'oswin'};F=sensc;for noliste=1:length(liste);hold on;t=liste{noliste};for i=.1:.1:.5;subplot(221);title(t);hold on;betiso=feval(F.iso.helperfunctions.baldev02isotherms,t,feval(F.iso.helperfunctions.baldev02STARCHandPE,t,i));betiso.awminawmax(1:2)=[.01 .9];betiso.awminawmax(3)=-1/30;F.iso.S=feval(F.iso.helperfunctions.approximate_isotherm,betiso);plot([-.1,F.iso.S.aws],[i,F.iso.S.Cs/F.iso.S.Cs(end),i]);subplot(222);hold on;plot([-.1,F.iso.S.aws],[i/10,F.iso.S.Cs/100,i/10]);subplot(223);hold on;densiteLDPE=922;pvsat27=3564.9138;RsurMeau=461.495;plot(F.iso.S.aws,[i/100,pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K(2:end)/densiteLDPE,i/100]);subplot(224);hold on;plot(F.iso.S.aws,[i,[pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K(2:end)/densiteLDPE]./max([pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K(2:end)/densiteLDPE]),i]);end;end;
%i;// clear;close all;liste={'bet' 'smith' 'halsey' 'caurie' 'oswin'};F=sensc;for noliste=1:length(liste);figure(noliste);hold on;t=liste{noliste};for i=.1:.1:.5;subplot(221);title(t);hold on;betiso=feval(F.iso.helperfunctions.baldev02isotherms,t,feval(F.iso.helperfunctions.baldev02STARCHandPE,t,i));betiso_awminawmax(1:2)=[.01 .9];betiso.awminawmax(3)=-1/30;F.iso.S=feval(F.iso.helperfunctions.approximate_isotherm,betiso);plot([-.1,F.iso.S.aws],[i,F.iso.S.Cs/F.iso.S.Cs(end),i]);subplot(222);hold on;plot([-.1,F.iso.S.aws],[i/10,F.iso.S.Cs/100,i/10]);subplot(223);hold on;densiteLDPE=922;pvsat27=3564.9138;RsurMeau=461.495;plot(F.iso.S.aws,[i/100,pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K(2:end)/densiteLDPE,i/100]);subplot(224);hold on;plot(F.iso.S.aws,[i,[pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K(2:end)/densiteLDPE]./max([pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K(2:end)/densiteLDPE]),i]);end;end;
%i;// clear;close all;liste={'bet' 'smith' 'halsey' 'caurie' 'oswin'};F=sensc;for noliste=1:length(liste);figure(noliste);hold on;t=liste{noliste};for i=.1:.1:.5;subplot(221);title(t);hold on;betiso=feval(F.iso.helperfunctions.baldev02isotherms,t,feval(F.iso.helperfunctions.baldev02STARCHandPE,t,i));betiso_awminawmax(1)=betiso.awminawmax(1)-.05;betiso.awminawmax(3)=-1/30;F.iso.S=feval(F.iso.helperfunctions.approximate_isotherm,betiso);plot([-.1,F.iso.S.aws],[i,F.iso.S.Cs/F.iso.S.Cs(end),i]);subplot(222);hold on;plot([-.1,F.iso.S.aws],[i/10,F.iso.S.Cs/100,i/10]);subplot(223);hold on;densiteLDPE=922;pvsat27=3564.9138;RsurMeau=461.495;plot(F.iso.S.aws,[i/1000,pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K(2:end)/densiteLDPE,i/1000]);subplot(224);hold on;plot(F.iso.S.aws,[i,[pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K(2:end)/densiteLDPE]./max([pvsat27/RsurMeau/(273.15+27)*100./F.iso.S.K(2:end)/densiteLDPE]),i]);end;end;
%i;// appel sous matlab : cd K:/recherCHE/matlab/senscrep,F=sensc,betiso=feval(F.iso.helperfunctions.baldev02isotherms,'smith',feval(F.iso.helperfunctions.baldev02STARCHandPE,'smith',.4));betiso.awminawmax(3)=-2-1/5;F.iso.S=feval(F.iso.helperfunctions.approximate_isotherm,betiso);F.C0sol=5;F.L=.1;[resode,respde,resanalytic]=sensc('analytic',F,true,true)
%i;// appel sous matlab : cd K:/recherCHE/matlab/senscrep/ ; [a,b,c]=sensc('anaLYTic',sensc,true,true) : OK, marche 20080610.
%i;// appel sous matlab : cd K:/recherCHE/matlab/senscrep/ ; F=sensc;F.C0sol=3;[a,b,c]=sensc('anaLYTic',F,true,true) : OK, marche 20080611
%i;// appel sous matlab : cd K:/recherCHE/matlab/senscrep/ ; F=sensc;F.C0sol=3;[a,b,c]=sensc('anaLYTic',F,true,true) : OK, marche 20080611
%i;// appel sous matlab : cd K:/recherCHE/matlab/senscrep/ ; F=sensc;F.C0sol=3;F.adimensionne.C0flu=.1;[a,b,c]=sensc('anaLYTic',F,true,true) : marche pas.
%i;// appel sous matlab : cd K:/recherCHE/matlab/senscrep;F=sensc;F.iso.S.K=[1 1];F.iso.S.Cs=[0 10];F.L=.1;F.Bi=1000;[resode,respde,resanalytic]=sensc('anaLYTIc',F,true,true) : marche 20080701, désavantage sagiv.

%  // sens pour sensibilité, et c pour corner. 
%  // maintenant, l'isotherme fluide est supposé etre c=PvsatsurRT*aw, et l'isotherme solide est supposé etre Cs(1)+sum(max([diff(Cs),eps^(-2)],K.*min(0,aw-aw0-[0,cumsum(diff(Cs)./K(1:end-1))]))), avec Cs, alpha, K des nombres réels vérifiant les contraintes diff(Cs)(i)>0, K(i)>0. TODO BUG tres mal gere si aw0 n'est pas nul.
%  // procedure faisant l'homogeneisation aux parametres non-homogeneises de l'utilisateur.
% on soustrait a la concentration dans le fluide la concentration initiale de ce fluide.
%  // errcatch(-1,'pause');exec('danieldebug.m');exec('sensc.m');F=sensc();sensc('analytic',F)


% TO BE MODIFIED ACCORDING TO THE PROCESSOR SPEED
global timeout_ode timeout_pde
%i;// for i in *.m; do echo '1,$s:endfunction:%&:UwUhUq' | tr U '\012' | ed $i; done
%i;// for i in *.m; do echo '1,$s:^\#:%&:UwUhUq' | tr U '\012' | ed $i; done
%i;// for i in *.m; do echo '1,$s:!=:~=:gUwUhUq' | tr U '\012' | ed $i; done
%i;// for i in *.m; do echo '1,$s:endif:end:gUwUhUq' | tr U '\012' | ed $i; done
%i;// for i in *.m; do echo '1,$s:!:~:gUwUhUq' | tr U '\012' | ed $i; done
%i;// for i in *.m; do echo '1,$s:  *(:(:gUwUhUq' | tr U '\012' | ed $i; done
%i;// for i in *.m; do echo '1,$s:usage:error:gUwUhUq' | tr U '\012' | ed $i; done
%i;// for i in *.m; do echo '1,$s:":'"'"':gUwUhUq' | tr U '\012' | ed $i; done
%i;// et ca continue plein de différences de syntaxe. 
%i;// addpath('octavestatpackage'); marche pas.
scilab=0;%i;scilab=1 // compatibilité scilab et matlab.
%i;// on n'a pas le package statistique de matlab.
if scilab==0
  if strcmp(version,'6.5.0.180913a (R13)');
    addpath('../scilabstatpackage'); % je préfere utiliser le package statistique de scilab, appelé directement
  end
end
timeout_ode		= 1e8; % s
timeout_pde		= 600e8; % s
%i;// [solode,solpde]=sensc('pde',u,1,1)
%i;// Definitions
isothermea1coin         = inline('iso.C_a+iso.K*min(aw,(iso.C_b-iso.C_a)/iso.K)+max(0,aw-(iso.C_b-iso.C_a)/iso.K)*iso.K*iso.alpha','aw','iso');
invisothermea1coin         = inline('(min(C,iso.C_b)-iso.C_a)/iso.K+max(0,C-iso.C_b)/iso.K/iso.alpha','C','iso'); % s = scale factor
GP          = inline('iso.PvsatsurRT*aw','aw','iso'); %i;// PvsatsurRT = scale factor, 4e-4 dans la realite, 1 apres eomogeneisation.
invGP          = inline('1/iso.PvsatsurRT*C','C','iso'); %i;// PvsatsurRT = scale factor, 4e-4 dans la realite, 1 apres homogeneisation.
%i;// xdefault	= mesh1D(50,[0 1],[.4 .2]);
xdefault=[0:.02:.28,.3:.03:.72,.75:.01:1]';
options		= odeset('RelTol',1e-5,'AbsTol',1e-6,'Stats','yes','Initialstep',1e-8,'Maxstep',.1,'Maxorder',4,'Events',@centiC0);
Fdefault	= 	struct(...
  'x',			xdefault,...
  'Bi',			10,...// Biot [hm.L1/D]
  'iso',		struct(...
    'S',		  struct('conc',isothermeacoin,'invconc',invisothermeacoin,'aw0',0,'Cs',[.1,.2],'K',[.5,1]),...
    'F',		  struct('conc',GP,'invconc',invGP,'PvsatsurRT',1)...
  ),...
  'L',			1/100,...// dimensionless thickness Lsol/Lflu [L1/L2]
  'adimensionne',	struct('L1',1,'Diff',1,'C0flu',0), ...// parametres qu'il faut adimensionner avant de faire tourner sensc.
  'C0sol',	1,...// dimensionless initial concentration [] , C0sol-C0flu, normalisé comme on veut ici.
  'xcrit',	0,...// dimensionless critical front position []
  'options',	options...
					);%i;// 
%i;// F.iso.S=struct(    'S',		  struct('conc',isothermea1coin,'invconc',invisothermea1coin,'C_a',.1,'C_b',.2,'alpha',2,'K',1));
Fdefault.t	= [0:.00001:.005 .01:.01:.1 .11:.1:5]';
Fdefault.t	= [0:.0005:.05 .06:0.01:.2 .22:0.02:.5 .55:.05:2 2.1:.1:5 6:1:10]';
Fdefault.iso.helperfunctions.baldev02STARCHandPE=@baldev02STARCHandPE;
Fdefault.iso.helperfunctions.baldev02isotherms=@baldev02isotherms;
Fdefault.iso.helperfunctions.approximate_isotherm=@approximate_isotherm;
method		= 'cubic';
Bi_critode	= 7e4;
Bi_maxode	= 1e18;
BiKminode	= 0.1;
MaxStepmax	= 10;
nstepchoice = 200;
n_default	= 1e4;
if Fdefault.Bi<10, Fdefault.t = Fdefault.t/(Fdefault.Bi/10); end
%i;// compatibilite avec vieux matlab
false=0;
true=1;
ploton_default = false;
dispon_default = false;
geometry_default = 0; %i;// slab


%i;[nargout,nargin]=argn(0);// Input control
if nargin<1 & ~nargout, methode = 'auto'; elseif nargin<1,	methode = 'init'; end
if nargin<2, F = Fdefault; end
%i;// adimensionnalisation. On commence par stocker les valeurs précédentes de ces paramètres.
verifieradimens=0; %i;// si ca vaut 1, ce qui ne doit pas nécessairement être adimensionalisé ne le sera pas. Utile pour vérifier certains paramètres.

F.adimensionne.nonadimensionne=F;
%i;// on enleve L1 et Diff.  % TODO verifier les lignes suivantes % theoreme vashy buckingham. les activités restent constantes.
F.t=F.t*F.adimensionne.Diff/F.adimensionne.L1^2;
first=1;%i;first=[]
for st={'F.C0sol', 'F.iso.S.Cs', 'F.adimensionne.C0flu', 'F.iso.S.K', 'F.iso.F.PvsatsurRT'}
  evalc(strcat([st{first},'=',st{first},'*F.adimensionne.L1;']));
end
F.x=F.x/F.adimensionne.L1;
%i;// ici, F.Bi contient le coefficient de transfert de masse (m/s)
F.Bi=F.Bi*F.adimensionne.L1/F.adimensionne.Diff ;
%i;// ici, F.Bi contient maintenant le nombre de biot.
F.adimensionne.L1=1;
F.adimensionne.Diff=1;
%i;// on enleve C0flu.
F.iso.S.aw0=F.iso.S.aw0-F.adimensionne.C0flu/F.iso.F.PvsatsurRT;
F.adimensionne.interne.C0flu=F.adimensionne.C0flu;
F.adimensionne.C0flu=0;
if verifieradimens==0
  %i;// on enleve Cs(1). On ne change pas les activités.
  F.C0sol=F.C0sol-F.iso.S.Cs(1);
  F.adimensionne.interne.Cs1=F.iso.S.Cs(1);
  F.iso.S.Cs=F.iso.S.Cs-F.iso.S.Cs(1);
  %i;// on enleve C0sol, qui est fixé à un. On ne change pas les activités.
  F.iso.S.Cs=F.iso.S.Cs/F.C0sol;
  F.iso.S.K=F.iso.S.K/F.C0sol;
  F.iso.F.PvsatsurRT=F.iso.F.PvsatsurRT/F.C0sol;
  F.adimensionne.interne.C0sol=F.C0sol;
  F.C0sol=1;
  %i;// on enleve PvsatsurRT, ce qui change les activités.
  F.iso.S.K=F.iso.S.K/F.iso.F.PvsatsurRT;
  F.adimensionne.interne.PvsatsurRT=F.iso.F.PvsatsurRT;
  F.iso.F.PvsatsurRT=1;
end
%i;// si l'adimensionnalisation a été demandée par l'utilisateur, on renvoie F avec les corrections faites.
if strcmp(methode,'adimensionaliselesparametres')
  resode=F;
  %i;// l'utilisateur peut voir les paramètres non adimensionnés dans resode.adimensionne.nonadimensionne, qui sera automatiquement écrasé au prochain lancement.
  return
end
if nargin<3, ploton = []; if ~nargout, ploton=true,end, end
  figreset=true;
if isa(ploton,'double')
  if mod(ploton,4)>1
    figreset=false
  end
  ploton=mod(ploton,2)==1
end
if nargin<4, dispon = []; end
if ~isfield(F,'autotime'), F.autotime = 'on'; end
if ~isfield(F,'n'), F.n = n_default; end
if ~isfield(F,'geometry'), F.geometry = geometry_default; end
if F.geometry>0 & (strcmp(methode,'ode') | strcmp(methode,'odepde'))
    if dispon, disp('... switch to PDE methode for non cartesian geometries'), end
    methode = 'pde';
end
if isempty(ploton), ploton=ploton_default; end
if isempty(dispon), dispon=dispon_default; end
doswitch=1;
while doswitch==1;
  doswitch=0;
  switch lower(methode);%i;select lower(methode);// compatibilité scilab matlab.
    case 'auto',	ploton = true; methode = 'all'; if dispon, disp('... auto methode'), end
	    %i;// DG20070509 pourquoi la ligne suivante ne finit-elle pas par "end, return, end" ? parce que les ends ci-dessus terminent des if ...
    case 'init',	resode = Fdefault; if dispon, disp('... init methode'), end, return
    case 'ode',		if dispon, disp('... ode methode'), end
    case 'pde',		if dispon, disp('... pde methode'), end
    case 'analytic',		if dispon, disp('... analytic methode'), end
    case 'odepde',	if nargout>1, if dispon, disp('... pde+ode methode'), end, else, methode = 'ode'; doswitch=1, end
    case 'all',	if nargout>2, if dispon, disp('... pde+ode+analytic methode'), end, else, methode = 'odepde'; doswitch=1, end
    otherwise;%i;else
      warning(strcat([methode,' and ',lower(methode),' are unknown command to sensc']))
      keyboard
      error('Sensc cannot repair input error');
  end
end
if ploton == true & figreset == true
  close all hidden
end
if strcmp(lower(F.autotime),'on')
	ti		= linspace(min(F.t),max(F.t),F.n)';
else
	ti		= F.t;
end
MaxStepmax = max(MaxStepmax,ceil(F.t(end)/nstepchoice));

%i;// update mesh to geometry
if F.geometry>0 & F.x(1)==0
    F.x = F.x + F.x(end)/F.L;
end

%i;// Solution
if any(strcmp(lower(methode),{'all' 'analytic'}))
    %i;// isothermeacoin         = inline('iso.C_a+iso.K*min(aw,(iso.C_b-iso.C_a)/iso.K)+max(0,aw-(iso.C_b-iso.C_a)/iso.K)*iso.K*iso.alpha','aw','iso');
    if ~inlinecmp(F.iso.S.conc,isothermeacoin)|~inlinecmp(F.iso.S.invconc,invisothermeacoin)
      error('l''isotherme solide doit avoir un coin et un seul');
    end
    if min(diff(F.iso.S.Cs))<0
      error('l''isotherme solide est mal paramétré, F.iso.S.Cs doit etre croissant.');
    end
    if ~inlinecmp(F.iso.F.conc,GP)|~inlinecmp(F.iso.F.invconc,invGP)
      error('l''isotherme liquide doit etre linéaire (henry)');
    end
    if F.iso.F.conc(F.iso.S.invconc(F.C0sol,F.iso.S),F.iso.F)<F.adimensionne.C0flu
      error('the analytic solution of sensc assumes that csbar will decrease');
    end
    dbstop if error
    %i;// ICI, on fait comme si rien n'avait ete adimensionne.
    %i;// coo premier isotherme.
    cs0=F.C0sol;
    cf0=F.iso.F.conc(F.iso.S.invconc(cs0,F.iso.S),F.iso.F);
    piso=max(find(F.iso.S.Cs<cs0));
    K0=F.iso.F.PvsatsurRT/F.iso.S.K(piso);
    %i;// 1+nblastterms, car il en faut au moins nblastterms pour calculer l'erreur dans Clast20 ...
    nblastterms=10; if isfield(F,'nblastterms'),nblastterms=F.nblastterms,end;
    moneps=1e-2;%i;// indication pour parametrer l'erreur globale en sortie.
    nmax=ceil(1+nblastterms+10*log(1+1/moneps)/sqrt(min(Fdefault.t(find(Fdefault.t>0)))));
    meth=[double(methode),zeros(1,10)+double('a')];
    if meth(4)<double('a'), if dispon, disp('Sagiv N=100'),end;
      nmax=100;
    end, if isfield(F,'nmax'),nmax=F.nmax;end;
    %i;// _t veut dire que la variable t parcourt un ensemble de valeurs (notation à la einstein) non nulles.
    %i;// __t veut dire fonction utilisée pour différentes valeurs de t.
    [lambda_n0,lambda_n0_cache,minpointfixep]=lambda_n__p(F.Bi,K0,F.L,nmax,meth(8)<double('a'));
    A_n00=A_n0__p(lambda_n0,lambda_n0_cache,F.Bi,K0,F.L);
    %i;// TODO, remplacer A_n0__p par A_n0__psurL et retarder le plus possible le produit par L.
    dbstop if error ;
    resanalytic=struct( ...// TODO ajouter l'energie mathematique ...
      'x',	F.x, ...
      'Cx',	zeros(1,size(F.x,1))*F.C0sol, ...// dans Cx, F.x', légende des colonnes, tC = légende des lignes
      'tC',	0, ...// le départ et les coins
      'fca4tC',	0, ...// juste avant passage des coins.
      'fcb4tC',	0, ...// flux sortant.
      'fc4tC',	0, ...// flux sortant
      'f4tC',	F.Bi*(cf0-F.adimensionne.C0flu), ...// flux sortant
      'Cm4tC',	F.C0sol, ...
      't',	F.t, ...
      'C',	zeros(size(F.t,1),1), ...
      'fc',	zeros(size(F.t,1),1), ...// fc :  perte totale.
      'fca',	zeros(size(F.t,1),1), ...// fca : perte à gauche, toujours nulle dans cette version de sensc.
      'fcb',	zeros(size(F.t,1),1), ...// fcb : perte à droite.
      'f',	zeros(size(F.t,1),1), ...// f ; flux à droite.
      'fD',	zeros(size(F.t,1),1) ...// flux divisé par Biot, je ne sais plus à quoi ça sert.
    );%i;// fait pour scilab pour contrer l'effet des ... ci-dessus sur les numéros de ligne. 
    %i;// pour t<0, pas de diffusion.
    segments=struct('a0p',F.C0sol,'a_np',0*A_n00,'lambda_np',lambda_n0,'lambda_np_cache',lambda_n0_cache,'minpointfixep',minpointfixep,'tp',-eps^(-2),'x',1,'csp',cs0,'cfp',cf0,'lastterms',nblastterms,'Kp',K0);
    p=2;
    segments(p)=struct('a0p',(K0*cs0+F.L*F.C0sol+F.adimensionne.C0flu-cf0)/(K0+F.L),'a_np',A_n00.*((cf0-F.adimensionne.C0flu-K0*cs0)/F.L+K0/F.L*F.C0sol),'lambda_np',lambda_n0,'lambda_np_cache',lambda_n0_cache,'minpointfixep',minpointfixep,'tp',0,'x',1,'csp',cs0,'cfp',cf0,'lastterms',nblastterms,'Kp',K0);
    debug=struct2cell(segments(2));if any(isnan(vertcat(debug{1:3,5:end})));error('nan detected for p=2');end;
    if meth(4)<double('a'), if dispon, disp('Sagiv init'),end;
      segment=segments(p);
      segments_energy=segments;
      Bn=(segment.csp-segment.a0p)*(sin(segment.lambda_np_cache.mod2pi))./segment.lambda_np;
      if meth(6)<double('a')
	[segments(p).a_np,segments(p).lambda_np_cache.normu]=analytic_inverse_Asagiv(segment,F.L,F.Bi,Bn);
	if meth(7)<double('a')
	  ttest=10.^(log10(eps):.1:0);
	  [conc,errconc]=sensca(ttest,1,segments(p));
	  [conc_energy,errconc_energy]=sensca(ttest,1,segments_energy(p));
	  resanalytic.err.sensca=struct('x',log10(ttest),'y',log10(eps^5+abs([conc,conc_energy]-F.C0sol)),'yerr',log10(eps^5+[abs(errconc(:,end))./exp(log(abs(eps^10+errconc(:,1))./abs(eps^10+errconc(:,end)))/(nblastterms-1)),abs(errconc_energy(:,end))./exp(log(abs(eps^10+errconc_energy(:,1))./abs(eps^10+errconc_energy(:,end)))/(nblastterms-1))]),'legend',struct('legend',{'sagiv ameliore','energy','sagiv errconv','energy errconv'}),'title','log10|C-1| at x=1 near log10(t=0)');
	  [conc,errconc]=senscam(ttest,segments(p));
	  [conc_energy,errconc_energy]=senscam(ttest,segments_energy(p));
	  resanalytic.err.senscam=struct('x',log10(ttest),'y',log10(eps^5+abs([conc,conc_energy]-F.C0sol)),'yerr',log10(eps^5+[abs(errconc(:,end))./exp(log(abs(eps^10+errconc(:,1))./abs(eps^10+errconc(:,end)))/(nblastterms-1)),abs(errconc_energy(:,end))./exp(log(abs(eps^10+errconc_energy(:,1))./abs(eps^10+errconc_energy(:,end)))/(nblastterms-1))]),'legend',struct('legend',{'sagiv ameliore','energy','sagiv errconv','energy errconv'}),'title','log10|Cm-1| near log10(t=0)')
	  [conc,errconc]=dsenscamsurdt(ttest,segments(p));
	  [conc_energy,errconc_energy]=dsenscamsurdt(ttest,segments_energy(p));
	  resanalytic.err.dsenscamsurdt=struct('x',log10(ttest),'y',log10(eps^5+abs([conc,conc_energy]+F.Bi*(cf0-F.adimensionne.C0flu))),'yerr',log10(eps^5+[abs(errconc(:,end))./exp(log(abs(eps^10+errconc(:,1))./abs(eps^10+errconc(:,end)))/(nblastterms-1)),abs(errconc_energy(:,end))./exp(log(abs(eps^10+errconc_energy(:,1))./abs(eps^10+errconc_energy(:,end)))/(nblastterms-1))]),'legend',struct('legend',{'sagiv ameliore','energy','sagiv errconv','energy errconv'}),'title','log10|f-f(t=0)| near log10(t=0)')
	  if ploton
	    figure(19);
	    subplot(231);
	    hold on
	    plot(resanalytic.err.sensca.x,resanalytic.err.sensca.y,resanalytic.err.sensca.x,resanalytic.err.sensca.yerr);
	    title(resanalytic.err.sensca.title);
	    legend(resanalytic.err.sensca.legend.legend);
	    subplot(233);
	    hold on
	    plot(resanalytic.err.senscam.x,resanalytic.err.senscam.y,resanalytic.err.senscam.x,resanalytic.err.senscam.yerr);
	    title(resanalytic.err.senscam.title);
	    legend(resanalytic.err.senscam.legend.legend);
	    subplot(234);
	    hold on
	    plot(resanalytic.err.dsenscamsurdt.x,resanalytic.err.dsenscamsurdt.y,resanalytic.err.dsenscamsurdt.x,resanalytic.err.dsenscamsurdt.yerr);
	    title(resanalytic.err.dsenscamsurdt.title);
	    legend(resanalytic.err.dsenscamsurdt.legend.legend);
	    if dispon
	      %i;// keyboard
	    end
	  end
	end
      else
        Asagiv=((eps^5+sin(kronenlogs(segment.lambda_np_cache.mod2pi,-segment.lambda_np_cache.mod2pi')))./(eps^5+kronenlogs(segment.lambda_np,-segment.lambda_np'))+(eps^5+sin(kronenlogs(segment.lambda_np_cache.mod2pi,segment.lambda_np_cache.mod2pi')))./(eps^5+kronenlogs(segment.lambda_np,segment.lambda_np')))/2;
	if meth(8)<double('a')
	  segments(p).a_np=inv(Asagiv)*Bn;%i;//  est TRES mauvais.
	else
	  segments(p).a_np=Asagiv\Bn;
	end
	if dispon
	  disp('conditionnement de Asagiv')
	  resanalytic.cond=cond(Asagiv);
	  resanalytic.rcond=rcond(Asagiv);
	  disp(resanalytic.cond);
	  disp(resanalytic.rcond);
	end
	if ploton
	  [s,v,d]=svd(Asagiv);
	  figure(19);
	  title('Comparison of Sagiv to energy')
	  subplot(231);
	  plot(diag(v));
	  title('spectrum of A')
	  subplot(233);
	  title('Test of residual inverse')
	  hold on
	  plot(Asagiv*segments(p).a_np)
	  plot(Bn,'g');
	  plot(log10(eps^3+abs(Asagiv*segments(p).a_np)),'c')
	  plot(log10(eps^3+abs(Bn)),'k');
	  plot(log10(eps^3+abs(Asagiv*segments(p).a_np-Bn)),'r');
	  legend('A*bn','Bn','log10 of A*bn','log10 of Bn','log10 of difference')
	end
      end
      if ploton
        figure(19);
	subplot(232);
        plot(log10(eps^5+abs([segments(p).a_np-segments_energy(p).a_np,segments(p).a_np,segments_energy(p).a_np])));
	legend('a_np error','a_np sagiv','a_np energy')
	title('log10 a_np and error')
      end
    end
    %i;// segment.x est l'endroit où le flux est important.
    %i;// si cs1 est entre cs0 et a00 :
    %i;// ATTENTION si la concentration en 1 remonte, pas de detection qu'on rechange d'isotherme !!
    while piso>1
      segment=segments(p);
      cspplus1=F.iso.S.Cs(piso);
      piso=piso-1;
      if (segment.a0p-cspplus1)*(cspplus1-segment.csp)>0
        cfpplus1=F.iso.F.conc(F.iso.S.invconc(cspplus1,F.iso.S),F.iso.F);
	Kpplus1=F.iso.F.PvsatsurRT/F.iso.S.K(piso);
	tplot=log(eps^.5+F.t);
	yplot=trouvecoin(tplot,segment,cspplus1);
	if ploton
	  figure(12)
	  clf
	  title('cs evolution without corner');
	  plot(tplot/log(10),cspplus1+0*F.t,'k')
	  hold on
	  plot(tplot/log(10),yplot+cspplus1)
	  legende=0;
	  legende.text='cs at corner';
	  legende(2).text='evolution';
	  if meth(5)<double('a')
	    plot(tplot/log(10),trouvecoin(tplot,segments_energy(p),cspplus1)+cspplus1,'g')
	    plot(tplot/log(10),log10(eps^5+abs(trouvecoin(tplot,segments_energy(p),cspplus1)-yplot)),'r')
	    legende(2).text='evolution with sagiv';
	    legende(3).text='evolution with energy';
	    legende(4).text='log10 evolution error';
	  end
	  legend(legende.text);
	end, if dispon, disp([strcat('Locating time at ',num2str(piso),'th corner (which is the ',num2str(p-1),'th corner we cross since beginning) ...')]),end; if dispon, DoDisplay='on'; else, DoDisplay='off'; end;
	if max(diff(yplot))>0
	  warning('Sensc gère mal l''adsorption, voir la condition if ci-dessus , plus ci-dessous cspplus1 devrait etre remplace par cs, couple avec for cs=[cspplus1,csp]. Verification supplémentaires en cours.');
	  if min(yplot(find(diff(yplot)>0)))<0 & max(yplot(find(diff(yplot)>0)))>0
	    error('this is at a crossing point')
	  end
	  if min(yplot(find(diff(yplot)>0)))<segment.csp-cspplus1 & max(yplot(find(diff(yplot)>0)))>segment.csp-cspplus1
	    error('l isotherme remonte un coin.')
	  end
	end
	%i;// avec l'arobas de matlab, ce serait plus propre, mais infaisable en scilab, et en vieux matlab.
	%i;// struct.cs1=0; [logt1,err]=lsqcurvefit(@trouvecoin,0,segment,cs1,[],[],optimset('Display',DoDisplay)); est une version plus lente et plus longue de la ligne suivante.
	%i;// partir de 0 plante lorsque c'est trop plat, genre lorsque biot est grand.
	[logt1,err]=fsolve(@trouvecoin,tplot(max([1;find(yplot(:)>0)])),optimset('Display',DoDisplay),segment,cspplus1);
	tpplus1=exp(logt1)+segment.tp;
	if abs(err)>eps^.1
	  warning('fsolve failed in sensc, or Bi is too high. Please enlarge F.nmax ...');
	  if tpplus1>eps^.4
	    warning('fsolve failed in sensc');
	  end
	end
	if dispon
	  disp([tpplus1 exp(logt1) cspplus1 segment.a0p segment.csp]) 
	end
	if meth(5)<double('a')
	  [logt1_energy,err]=fsolve(@trouvecoin,0,optimset('Display',DoDisplay),segments_energy(p),cspplus1);
	  tpplus1_energy=exp(logt1_energy)+segments_energy(p).tp;
	  if dispon
	    disp('time at corner:[error at log in step p,cumulated error,energy time,sagiv time]')
	    disp(mat2str([logt1_energy-logt1,tpplus1_energy-tpplus1,tpplus1_energy,tpplus1]));
	  end
	end
	resanalytic.tC(end+1,1)=tpplus1;
	resanalytic.Cx(end+1,:)=sensca(tpplus1,resanalytic.x,segment)';
	resanalytic.fca4tC(end+1,1)=0;
	resanalytic.fcb4tC(end+1,1)=-primitivetgradientxbordanalytic(tpplus1,segments);
	%i;// il ne fallait pas multiplier ceci par F.Bi !!!!
	resanalytic.fc4tC(end+1,1)=F.C0sol-senscam(tpplus1,segment);
	resanalytic.f4tC(end+1,1)=F.Bi*(cfpplus1-F.adimensionne.C0flu-F.L*resanalytic.fc4tC(end,1));
	resanalytic.Cm4tC(end+1,1)=senscam(tpplus1,segment);
	[lambda_npplus1,lambda_npplus1_cache,minpointfixep]=lambda_n__p(F.Bi,Kpplus1,F.L,nmax,meth(8)<double('a'));
	a_npplus1=zeros(deal(size(lambda_npplus1)));
        if meth(5)<double('a') | meth(4)>=double('a')
	  tpplus1_sagiv=tpplus1;
	  if meth(5)<double('a')
	    tpplus1=tpplus1_energy;
	    segment=segments_energy(p);
	  end
	  A_n0pplus1=A_n0__p(lambda_npplus1,lambda_npplus1_cache,F.Bi,Kpplus1,F.L);
	  %i;// TODO, remplacer A_n0__p par A_n0__psurL et retarder le plus possible le produit par L.
	  if meth(1)<double('a')
	    %i;// b_n1+A_n0pplus1*a00
	    a_npplus1=a_npplus1+A_n0pplus1*(cfpplus1-F.adimensionne.C0flu-F.L*F.C0sol)/F.L+A_n0pplus1*segment.a0p;
	  end
	  if meth(2)>=double('a')
	    %i;// b_n1prime+A_n01prime*a00
	    a_npplus1=a_npplus1+A_n0pplus1*(cfpplus1-F.adimensionne.C0flu-F.L*F.C0sol-Kpplus1*cspplus1)/F.L+A_n0pplus1*(F.L+Kpplus1)/F.L*segment.a0p;
	  end
	  for n=1:nmax;
	    An_mpplus1=zeros(nmax,1);
	    if meth(1)<double('a')
	      An_mpplus1=An_mpplus1+A_n0pplus1(n).*(Kpplus1*segment.lambda_np.*exp(-(tpplus1-segment.tp)*segment.lambda_np.^2))./(2*F.L*sin(lambda_npplus1_cache.mod2pi(n))).*(2*F.L*sin(segment.lambda_np_cache.mod2pi)*sin(lambda_npplus1_cache.mod2pi(n))./(Kpplus1*segment.lambda_np.^2)+(eps^5+sin(segment.lambda_np_cache.mod2pi-lambda_npplus1_cache.mod2pi(n)))./(eps^5+(segment.lambda_np-lambda_npplus1(n)))-(eps^5+sin(segment.lambda_np_cache.mod2pi+lambda_npplus1_cache.mod2pi(n)))./(eps^5+(segment.lambda_np+lambda_npplus1(n))));
	    end
	    if meth(2)>=double('a')
	      %i;// TODO developper les sinus et utiliser les caches avec valeur de cosinus et de sinus.
	      An_mpplus1=An_mpplus1+A_n0pplus1(n).*(Kpplus1*lambda_npplus1(n).*exp(-(tpplus1-segment.tp)*segment.lambda_np.^2))/(2*F.L*sin(lambda_npplus1_cache.mod2pi(n))).*(2*F.L*sin(segment.lambda_np_cache.mod2pi)*sin(lambda_npplus1_cache.mod2pi(n))./(Kpplus1*segment.lambda_np*lambda_npplus1(n))+(eps^5+sin(segment.lambda_np_cache.mod2pi-lambda_npplus1_cache.mod2pi(n)))./(eps^5+(segment.lambda_np-lambda_npplus1(n)))+(eps^5+sin(segment.lambda_np_cache.mod2pi+lambda_npplus1_cache.mod2pi(n)))./(eps^5+(segment.lambda_np+lambda_npplus1(n))));
	    end
	    %i;// le produit matriciel pour changer de base.
	    a_npplus1(n)=a_npplus1(n)+An_mpplus1'*segment.a_np;
	  end
	  if meth(1:2)==double('An')
	    a_npplus1=a_npplus1/2;
	  end
	  if meth(1:2)==double('aN')
	    error('il faut au moins l une des deux méthodes pour le produit matriciel')
	  end
	  tpplus1=tpplus1_sagiv;
	end
	p=p+1;
	segments(p)=struct('a0p',(Kpplus1*cspplus1+F.L*F.C0sol+F.adimensionne.C0flu-cfpplus1)/(Kpplus1+F.L),'a_np',a_npplus1,'lambda_np',lambda_npplus1,'lambda_np_cache',lambda_npplus1_cache,'minpointfixep',minpointfixep,'tp',tpplus1,'x',1,'csp',cspplus1,'cfp',cfpplus1,'lastterms',nblastterms,'Kp',Kpplus1);
	debug=struct2cell(segments(2));if any(isnan(vertcat(debug{1:3,5:end})));error('nan detected for p=2');end;
        if meth(4)<double('a'), if dispon, disp('Sagiv corner'),end;
          segment=segments(p);
          Bsagiv=((eps^5+sin(kronenlogs(segment.lambda_np_cache.mod2pi,-segments(p-1).lambda_np_cache.mod2pi')))./(eps^5+kronenlogs(segment.lambda_np,-segments(p-1).lambda_np'))+(eps^5+sin(kronenlogs(segment.lambda_np_cache.mod2pi,segments(p-1).lambda_np_cache.mod2pi')))./(eps^5+kronenlogs(segment.lambda_np,segments(p-1).lambda_np')))/2;
	  if meth(5)<double('a')
	    segments_energy(p)=segments(p);
	    segments_energy(p).tp=tpplus1_energy;
	  end
	  Bn=Bsagiv*(exp(-(segment.tp-segments(p-1).tp)*segments(p-1).lambda_np.^2).*segments(p-1).a_np)+(segments(p-1).a0p-segment.a0p)*(sin(segment.lambda_np_cache.mod2pi))./segment.lambda_np;
	  if meth(6)<double('a')
	    [segments(p).a_np,segments(p).lambda_np_cache.normu]=analytic_inverse_Asagiv(segment,F.L,F.Bi,Bn);
	  else
	    Asagiv=((eps^5+sin(kronenlogs(segment.lambda_np_cache.mod2pi,-segment.lambda_np_cache.mod2pi')))./(eps^5+kronenlogs(segment.lambda_np,-segment.lambda_np'))+(eps^5+sin(kronenlogs(segment.lambda_np_cache.mod2pi,segment.lambda_np_cache.mod2pi')))./(eps^5+kronenlogs(segment.lambda_np,segment.lambda_np')))/2;
	    if meth(8)<double('a')
	      segments(p).a_np=inv(Asagiv)*Bn;%i;//  est TRES mauvais.
	    else
	      segments(p).a_np=Asagiv\Bn;
	    end
	    if ploton
	      [s,v,d]=svd(Asagiv);
	      figure(19);
	      subplot(234);
	      plot(diag(v));
	      title(strcat(['spectrum of Asagiv at ',num2str(p-1),'th corner']))
	    end
	  end
	  if ploton
	    if meth(5)<double('a')
              figure(19);
	      subplot(2,(p+2)*(p+3)/2,(p+2)*(p+3)/2+((p+3)*(p-2)/2+1:(p+2)*(p-1)/2));
	      a_npsagivfromenergy=segments_energy(p).a_np.*exp(-(segment.tp-segments_energy(p).tp)*segments_energy(p).lambda_np.^2);
              plot(log10(eps^5+abs([segments(p).a_np-a_npsagivfromenergy,segments(p).a_np,a_npsagivfromenergy])));
	      legend('a_np error','a_np sagiv','a_np energy')
	      title(strcat(['log10 a_n and error after ',num2str(p-1),'th corner']))
	    end
	  end
	  debug=struct2cell(segments(p));if any(isnan(vertcat(debug{1:3,5:end})));error('nan detected for p');end;
        end
      else
	piso=-1;
      end
    end
    scilab=0;%i;scilab=1
    resanalytic.functions.value_t_x_param=@sensca;
    resanalytic.functions.mean_timet_param=@senscam;
    resanalytic.functions.moinsflux_timet_param=@dsenscamsurdt;
    if ~strcmp(version,'6.5.0.180913a (R13)') & scilab==0
      evalc('resanalytic.functions.value_t_x=@(t,x)sensca(t,x,segments);');
      evalc('resanalytic.functions.mean_timet=@(t)senscam(t,segments);');
      evalc('resanalytic.functions.flux_timet=@(t)dsenscamsurdt(t,segments);');
    end
    tpplus1=eps^(-.5);
    resanalytic.tC(end+1,1)=tpplus1;
    resanalytic.Cx(end+1,:)=sensca(tpplus1,resanalytic.x,segments)';
    resanalytic.fca4tC(end+1,1)=0;
    resanalytic.fcb4tC(end+1,1)=-primitivetgradientxbordanalytic(tpplus1,segments);
    %i;// il ne fallait pas multiplier ceci par F.Bi !!!!
    resanalytic.fc4tC(end+1,1)=F.C0sol-senscam(tpplus1,segments);
    resanalytic.f4tC(end+1,1)=-dsenscamsurdt(tpplus1,segments);
    resanalytic.Cm4tC(end+1,1)=senscam(tpplus1,segments);
    for not=1:size(F.t,1);
      t=F.t(not);
      if t==0
        resanalytic.C(not)=F.C0sol;
	resanalytic.fc(not)=0;
	resanalytic.fca(not)=0;
	resanalytic.fcb(not)=0;
	resanalytic.f(not)=F.Bi*(cf0-F.adimensionne.C0flu);
      else
	resanalytic.C(not)=senscam(t,segments);
	resanalytic.fc(not)=F.C0sol-senscam(t,segments);
	resanalytic.fca(not)=0;
	resanalytic.fcb(not)=-primitivetgradientxbordanalytic(t,segments);
	%i;// il ne fallait pas multiplier fcb(not) par Bi !!!
	resanalytic.f(not)=-dsenscamsurdt(t,segments);
      end
    end
    if meth(3)<double('a') | ploton==true
      departetdautre=[];
      if meth(7)>double('a')
	departetdautre=kron(10.^(0:-2:log10(eps)*.5),[-1 1]);
	departetdautre=departetdautre(:);
	departetdautre=unique(sort([kronenlogs([segments(2:end).tp],departetdautre),kron([segments(2:end).tp],1+0*departetdautre)+kron(0*[segments(2:end).tp],departetdautre)]));
	departetdautre=departetdautre(find(departetdautre>0));
      end
      listet=[resanalytic.t;departetdautre];
      for not=1:size(listet,1);
	t=listet(not);
        resanalytic.Cx(end+1,:)=sensca(t,F.x,segments)';
      end
      %i;// on se debarasse des t doublons.
      [resanalytic.tC,ordre,choixdanst]=unique([resanalytic.tC;listet]);
      resanalytic.Cx=resanalytic.Cx(ordre,:);
      ind=substruct('()',{ordre});
      resanalytic.fca4tC=subsref([resanalytic.fca4tC;resanalytic.fca;0*departetdautre],ind);
      resanalytic.fcb4tC=subsref([resanalytic.fcb4tC;resanalytic.fcb;-primitivetgradientxbordanalytic(departetdautre,segments)],ind);
      resanalytic.fc4tC=subsref([resanalytic.fc4tC;resanalytic.fc;F.C0sol-senscam(departetdautre,segments)],ind);
      resanalytic.f4tC=subsref([resanalytic.f4tC;resanalytic.f;-dsenscamsurdt(departetdautre,segments)],ind);
      resanalytic.Cm4tC=subsref([resanalytic.Cm4tC;resanalytic.C;senscam(departetdautre,segments)],ind);
      %i;// ces quatres lignes sont équivalents :
      %i;// a=[a;b];a=a(ordre,:);
      %i;// [bidon,ind]=sort(ordre,:); a(ind)=[a;b];
      %i;// ind=1; ind(ordre)=1:length(ordre); a(ind)=[a;b]; attention pas de possibilité d'ajouter une elimination de duplicatas ici.
      %i;// a=subsref([a;b],substruct('()',{ordre,':'}));
    end
    resanalytic.parametres=segments;
    resanalytic.parametres_energy=[];
    if meth(5)<double('a')
      resanalytic.parametres_energy=segments_energy;
    end
    resanalytic.fD=resanalytic.f/F.Bi;
    if any(strcmp(lower(methode),{'analytic'}))
      resode=resanalytic;
      resode.Cx=[];
      resode.parametres=[];
      resode.parametres_energy=[];
      resode.functions=[];
      respde=resode;
    end
end
if any(strcmp(lower(methode),{'all' 'odepde' 'pde'}))
	if dispon, disp('SENSpde: solving...'), end
    %i;//  	if F.Bi*F.K>50, F.options.MaxOrder= 2; end
    errpde = 0;
    try
        tic, [sol,Fca,Fcb,tC,Cx,FcatC,FcbtC] = pdesolve1D(F.geometry,@pdediff_PDE,@pdediff_IC,@pdediff_BC,F.x,F.t,F.options,F);
    catch
        if findstr(lasterr,'erron')
            if dispon, disp('... TIME OUT'), end
            timeout_pdeold = timeout_pde; timeout_pde = 15;
            F.options = odeset(F.options,'bdf','on','maxorder',2,'maxstep',1e5);
            try
                tic, [sol,Fca,Fcb,tC,Cx,FcatC,FcbtC] = pdesolve1D(F.geometry,@pdediff_PDE,@pdediff_IC,@pdediff_BC,F.x,F.t,F.options,F);
            catch
                if findstr(lasterr,'erron'), if dispon, disp('... TIME OUT (2)'), end, else, if dispon, disp(sprintf('...unexpected error ''%s'' for F =',lasterr)), disp(F), end; end
                errpde = 1;
            end
            timeout_pde = timeout_pdeold;
        else
            if dispon, disp('... decrease maxstep'), end
            timeout_pdeold = timeout_pde; timeout_pde = 30;
            F.options.MaxStep = F.options.MaxStep/100;
            try
                tic, [sol,Fca,Fcb,tC,Cx,FcatC,FcbtC] = pdesolve1D(F.geometry,@pdediff_PDE,@pdediff_IC,@pdediff_BC,F.x,F.t,F.options,F);
            catch
                if findstr(lasterr,'erron'), if dispon, disp('... TIME OUT'), end, else, if dispon, disp(sprintf('...unexpected error ''%s'' for F =',lasterr)), disp(F); end, end
                errpde = 1;
            end
            timeout_pde = timeout_pdeold;        
        end
    end
    if ~errpde
        if dispon, disp(['SENSpde: end in ' num2str(toc) ' s']), end
        % Extraction
	C		= Cx(:,:,1);
        C1		= C(:,end);
        F.x     = F.x(:)';
        if F.geometry>0
            Cm		= trapz(F.x,C.*repmat(F.x,size(C,1),1).^F.geometry,2) / trapz(F.x,F.x.^F.geometry);
            fc		= pdeloss(F.x,C,0,F.geometry);
        else
            Cm		= trapz(F.x,C,2);
            fc		= pdeloss(F.x,C);
        end
        Cmi		= interp1(tC,Cm,ti,method);
        fci		= interp1(tC,fc,ti,method);
        C1i		= interp1(tC,C1,ti,method);
	%i;// une ligne de DG
        fi		= interp1((tC(1:end-1)+tC(2:end))/2,(fc(2:end)-fc(1:end-1)).*(tC(2:end)-tC(1:end-1)),ti,method);
        % fi		= F.Bi*F.K*C1i-F.Bi*F.L*fci;
        if fi(end)<0, fi = fi-min(fi); end
        %i;// 	fi(1)	= F.K*F.Bi*F.C0sol;
            respde.x = F.x;
            respde.Cx = C;
	    %i;// x = légende des colonnes, tC = légende des lignes
            respde.tC = tC;
	    %i;// l'erreur ci-dessous est basse frequence
            respde.fca4tC = FcatC;
            respde.fcb4tC = FcbtC;
	    %i;// l'erreur ci-dessous est haute frequence
            respde.fc4tC = fc;
            respde.Cm4tC = Cm;
            respde.t = ti;
            respde.C = Cmi;
            respde.fc = fci;
            respde.fca = Fca;
            respde.fcb = Fcb;
            respde.f = fi;
            respde.fD = fi/F.Bi;
    else
        if dispon, disp(['SENSpde is TIMEOUT => switch to SENSode']), end
    end
end
%i;// if any(strcmp(lower(methode),{'all' 'odepde' 'ode'})) marche plus car on a enlevé .K
if any(strcmp(lower(methode),{'odepde' 'ode'}))
    tic, if dispon, disp('SENSode: solving...'), end
    nsteptypical = ceil(F.t(end)/F.options.MaxStep);
    if (nsteptypical>nstepchoice) | (F.Bi*F.K<BiKminode)
        if dispon, disp(sprintf('... increase maxstep up to %0.2g',MaxStepmax)), end
        factor = ceil(MaxStepmax/F.options.MaxStep);
        F.options = odeset(F.options,'maxstep',MaxStepmax,'Initialstep',factor*F.options.InitialStep);
    elseif F.Bi*F.K<10*BiKminode
        step = MaxStepmax*BiKminode/(F.Bi*F.K);
        if dispon, disp(sprintf('... increase maxstep up to %0.2g',step)), end
        factor = ceil(step/F.options.MaxStep);
        F.options = odeset(F.options,'maxstep',step,'Initialstep',factor*F.options.InitialStep);
    end
    if F.Bi<Bi_maxode
        try
            if F.Bi>=Bi_critode
                F.options = odeset(F.options,'bdf','on','maxorder',2,'maxstep',.1);
                if dispon, disp('... bdf is ''on'''), end
            end
            tic
            [tC,Code]	= ode15s(@odediff,F.t,[F.C0sol],F.options,F);
        catch
            if findstr(lasterr,'erron')
                if dispon, disp('... TIME OUT'), end
                F.options = odeset(F.options,'maxstep',1e5,'Initialstep',.1);
                timeout_odeold = timeout_ode; timeout_ode = 10;
                tic
                try
                    [tC,Code]	= ode15s(@odediff,F.t,[F.C0sol],F.options,F);
                catch
                    if dispon, disp(sprintf('...unexpected error ''%s'' for F =',lasterr)), disp(F), end
                    if dispon, disp(sprintf('==>   STOP at %s', datestr(now))), end
                    tC		= F.t;
                    Code	= ones(size(tC));
                end
                timeout_ode = timeout_odeold;
            else
                if dispon, disp('... decrease maxstep'), end
                F.options.MaxStep = F.options.MaxStep/100;
                tic
                [tC,Code]	= ode15s(@odediff,F.t,[F.C0sol],F.options,F);
            end
        end
    else
        tC = F.t;
        Code = F.C0sol * ones(size(F.t));
    end
    if dispon, disp(['SENSode: end in ' num2str(toc) ' s']), end
    Codei	= interp1(tC,Code,ti,method);
    %fodei	= - ndf(ti,Codei);
    timeout_odeold = timeout_ode; timeout_ode = 10; tic
    fodei	= - odediff(Codei,F);
    timeout_ode = timeout_odeold;
end

% Outputs
avecresode=0;
avecrespde=0;
if any(strcmp(lower(methode),{'all' 'odepde' 'ode'}))
	if nargout>0
		resode.t = ti;
		resode.C = Codei;
		resode.f = fodei;
		resode.fD = fodei/F.Bi;
	end
end
if any(strcmp(lower(methode),{'all' 'odepde' 'pde'}))
  if nargout>0 & strcmp(lower(methode),'pde'), resode	= []; end
  if nargout>1
    if ~errpde
    else
      switch lower(methode);%i;select lower(methode);// compatibilité scilab matlab.
	case 'pde', error('resolution ratee, et resode n est pas encore au point pour ../sensc');%i;//respde = sensc('pde',F,ploton);
	case 'odepde', respde = resode;
      end
    end
  end
end
% desadimensionnalisation.
if verifieradimens==0
  % on remet PvsatsurRT
  F.iso.F.PvsatsurRT=F.adimensionne.interne.PvsatsurRT;
  F.iso.S.K=F.iso.S.K*F.iso.F.PvsatsurRT;
  % on remet C0sol.
  F.C0sol=F.adimensionne.interne.C0sol;
  for st=allongelistest({'resanalytic.Cx','resanalytic.fca4tC','resanalytic.fcb4tC','resanalytic.fc4tC','resanalytic.f4tC','resanalytic.Cm4tC','resanalytic.C','resanalytic.fc','resanalytic.fca','resanalytic.fcb','resanalytic.f','resanalytic.fD','resode.C', 'resode.f', 'resode.fD', 'respde.Cx', 'respde.fc4tC', 'respde.Cm4tC', 'respde.C', 'respde.fc', 'respde.fca', 'respde.fcb', 'respde.f', 'respde.fD','F.iso.S.Cs','F.iso.F.PvsatsurRT','F.iso.S.K'},{'a0p','a_np','csp','Kp'},resanalytic);
    evalc(strcat([st{first},'=',st{first},'*F.C0sol;']),'''pas grave'';');
  end
  %i;//    if meth(5)<double('a') parametres_energy reste adimensionalise !!!
  % on remet Cs(1).
  F.iso.S.Cs=F.iso.S.Cs+F.adimensionne.interne.Cs1;
  for st=allongelistest({'resanalytic.Cx','resanalytic.Cm4tC','resanalytic.C','resode.C','respde.Cx', 'respde.Cm4tC', 'respde.C','F.C0sol'},{'csp', 'a0p'},resanalytic);
    evalc(strcat([st{first},'=',st{first},'+F.iso.S.Cs(1);']),'''pas grave'';');
  end
end
% on remet C0flu.
F.adimensionne.C0flu=F.adimensionne.interne.C0flu;
F.iso.S.aw0=F.iso.S.aw0+F.adimensionne.C0flu/F.iso.F.PvsatsurRT;
for st=allongelistest({},{'cfp'},resanalytic);
  evalc(strcat([st{first},'=',st{first},'+F.adimensionne.C0flu;']),'''pas grave'';');
end
% on remet L1 et Diff.
F.adimensionne.L1=F.adimensionne.nonadimensionne.adimensionne.L1;
F.adimensionne.Diff=F.adimensionne.nonadimensionne.adimensionne.Diff;
for st=allongelistest({'resanalytic.tC','resanalytic.t','F.t','resode.t','respde.t','respde.tC'},{'tp'},resanalytic)
  evalc(strcat([st{first},'=',st{first},'/F.adimensionne.Diff*F.adimensionne.L1^2;']),'''pas grave'';');
end
for st=allongelistest({'resanalytic.x','F.x','respde.x'},{'x'},resanalytic)
  evalc(strcat([st{first},'=',st{first},'*F.adimensionne.L1;']),'''pas grave'';');
end
for st=allongelistest({},{'lambda_np'},resanalytic)
  evalc(strcat([st{first},'=',st{first},'/F.adimensionne.L1;']),'''pas grave'';');
end
for st=allongelistest({'resanalytic.Cx','resanalytic.Cm4tC','resanalytic.C','resode.C','respde.C','respde.Cm4tC','respde.Cx','F.C0sol', 'F.iso.S.Cs', 'F.adimensionne.C0flu', 'F.iso.S.K', 'F.iso.F.PvsatsurRT'},{'a0p','a_np','csp','Kp'},resanalytic)
  evalc(strcat([st{first},'=',st{first},'/F.adimensionne.L1;']),'''pas grave'';');
end
for st={'resanalytic.fD','resode.fD','respde.fD'}
  evalc(strcat([st{first},'=',st{first},'/F.adimensionne.L1^2;']),'''pas grave'';');
end
for st={'resanalytic.f4tC','resanalytic.f','resode.f','respde.f'}
  evalc(strcat([st{first},'=',st{first},'*F.adimensionne.Diff/F.adimensionne.L1^3;']),'''pas grave'';');
end
%i;// ici, F.Bi contient encore le nombre de biot.
F.Bi=F.Bi/F.adimensionne.L1*F.adimensionne.Diff;
%i;// ici, F.Bi contient maintenant le coeffcient de transfert de masse (m/s)
F.x=F.x*F.adimensionne.L1;
	
Fostring='t'; if F.adimensionne.Diff==1 & F.adimensionne.L1==1, Fostring='Fo'; end
Bistring='D'; if F.adimensionne.Diff==1 & F.adimensionne.L1==1, Bistring='Bi'; end
if ploton
	if any(strcmp(lower(methode),{'all' 'odepde' 'pde'}))
		figure(1)
		clf
		subplot(231), hold on, plot(F.x,respde.Cx), xlabel('x'), ylabel('C')
		subplot(232), hold on, plot(ti,Cmi,'b-'), xlabel(Fostring)
		subplot(233), hold on, plot(sqrt(ti),Cmi,'b-'), xlabel([Fostring,'^{1/2}'])
		subplot(234), hold on, plot(ti,fi,'b-'), ylabel('j')
		subplot(236), hold on, plot(Cmi,fi,'b-'), xlabel('C'), ylabel('j'), title('kpd')
		figure(2)
		clf
		subplot(121), hold on, plot(ti/F.Bi,Cmi,'k:'), xlabel([Fostring,'/Bi']), ylabel('C')
		subplot(122), hold on, plot(Cmi,1/F.Bi*fi,'k:'), xlabel('C'), ylabel('j/Bi')
		figure(3)
		clf
		maxnot=size(respde.tC,1),
		maxnox=size(F.x,1),
		tcolor=interp1(0:11,jetcolormap12(0),11*(0:maxnot-1)/(maxnot-.9),'linear');
		xcolor=interp1(0:11,jetcolormap12(0),11*(0:maxnox-1)/(maxnox-.9),'linear');
		subplot(224),hold on,for not=1:maxnot;plot(respde.Cx(not,:),F.x,'color',tcolor(not,:)),end,ylabel('x'),xlabel('C')
		subplot(221),hold on,for nox=1:maxnox;plot(respde.tC,respde.Cx(:,nox),'color',xcolor(nox,:)),end,xlabel(Fostring),ylabel('C')
		%i;// ensia ensbana initia agro : veronique a torena ?
		Cxpourdt=[respde.Cx,respde.tC.^2/2];
		dCxsurdt=spdiags((respde.tC(3:end)-respde.tC(1:end-2)).^(-1),0,maxnot-2,maxnot-2)*(Cxpourdt(3:end,:)-Cxpourdt(1:end-2,:));
		%i;// pour que interp1 n'est pas deux fois la meme abscisse.
		trouvezeros=find(diff(dCxsurdt(:,end))==0);
		dCxsurdt(trouvezeros+1,end)=min(dCxsurdt(trouvezeros,end)+eps*max(1,abs(dCxsurdt(trouvezeros,end))),(dCxsurdt(trouvezeros+2,end)+dCxsurdt(trouvezeros,end))/2);
		dCxsurdt=interp1(dCxsurdt(:,end),dCxsurdt(:,1:end-1),respde.tC(2:end-1),'linear');
		Cxpourdxdx=[respde.Cx;F.x'.^3/6];
		dCxsurdxgauche=(Cxpourdxdx(:,2:end)-Cxpourdxdx(:,1:end-1))*spdiags((F.x(2:end)-F.x(1:end-1)).^(-1),0,maxnox-1,maxnox-1); %i;// dCxsurdx decentre a droite
		ddCxsurdxdx=(dCxsurdxgauche(:,2:end)-dCxsurdxgauche(:,1:end-1))*spdiags((F.x(3:end)-F.x(1:end-2)).^(-1),0,maxnox-2,maxnox-2)*2; %i;// et on redérive. Ce n'est pas tellement centré :-(
		ddCxsurdxdx=interp1(ddCxsurdxdx(end,:),ddCxsurdxdx(1:end-1,:)',F.x(2:end-1),'linear')';
		errrelative=(ddCxsurdxdx(2:end-1,:)*F.adimensionne.Diff-dCxsurdt(:,2:end-1))./(eps^.5+max(abs(ddCxsurdxdx(2:end-1,:))*F.adimensionne.Diff,abs(dCxsurdt(:,2:end-1))));
		listetrait=kronenlogs([.2 .3 .5 .7],0*[-6:0]).*(10.^kronenlogs(0*[.2 .3 .5 .7],[-6:0]));
		subplot(223),hold on,contour(respde.tC(2:end-1)',F.x(2:end-1),errrelative',[-listetrait(end:-1:1),0,listetrait]);[h,c]=contour(respde.tC(2:end-1)',F.x(2:end-1),errrelative'*1e6,[-10.^[6:-1:0],0,10.^[0:6]]);clabel(h,c);xlabel(Fostring),ylabel('x')
		for not=1:maxnot;plot(respde.tC(not),max(F.x),'color',tcolor(not,:));end
		for nox=1:maxnox;plot(max(respde.tC),F.x(nox),'color',xcolor(nox,:));end
		fluxbord=(2*F.x(end)-F.x(end-1)-F.x(end-2))'*dCxsurdxgauche(1:end-1,end)+(F.x(end)-F.x(end-1))'*dCxsurdxgauche(1:end-1,end-1);
		fluxbord=-F.adimensionne.Diff*fluxbord./(2*F.x(end)-F.x(end-1)-F.x(end-2)+F.x(end)-F.x(end-1));
		fluxbordbis=(respde.fcb4tC(3:end)-respde.fcb4tC(1:end-2))./(respde.tC(3:end)-respde.tC(1:end-2));
		fluxbordter=(respde.fc(3:end)-respde.fc(1:end-2))./(respde.t(3:end)-respde.t(1:end-2));
		subplot(443),hold on,plot(respde.tC,fluxbord,respde.tC(2:end-1),fluxbordbis,respde.t(2:end-1),fluxbordter,respde.tC(2:end-1),fluxbord(2:end-1)-fluxbordbis,respde.t(2:end-1),interp1(respde.tC,fluxbord,respde.t(2:end-1),'linear')-fluxbordter),legend('flux par Cx','flux par -fcb4tC','flux par fc','difference Cx a -fcb4tC','difference Cx a fc'),xlabel(Fostring),ylabel('flux')
		danscolonne=1;%i;danscolonne='r';// pour scilab
		subplot(444),hold on,plot(respde.tC,respde.Cm4tC,respde.t,respde.C,respde.t,F.C0sol-respde.fc/F.adimensionne.L1,respde.t,interp1(respde.tC,respde.Cm4tC,respde.t,'linear')-F.C0sol+respde.fc/F.adimensionne.L1,respde.t,respde.C-F.C0sol+respde.fc/F.adimensionne.L1),legend('Cmoyen par Cm4tC','Cmoyen par C','Cmoyen par fc','difference Cm4tC a fc','difference C a fc'),xlabel(Fostring),ylabel('concentration')
		concentrationinterfaceliquide=F.adimensionne.C0flu-respde.Cm4tC*F.L+F.C0sol*F.L+fluxbord/F.Bi;
		Cs=[F.iso.S.Cs,F.C0sol]; if Cs(end)<Cs(end-1); Cs=Cs(1:end-1); end;
		subplot(447),hold on,plot(F.iso.F.invconc(concentrationinterfaceliquide,F.iso.F),respde.Cx(:,end),F.iso.S.invconc(Cs,F.iso.S),Cs,F.iso.S.invconc(Cs,F.iso.S),Cs,'*',F.iso.F.invconc(F.adimensionne.C0flu,F.iso.F),F.C0sol,'+',F.iso.S.invconc(F.Cx(end,end),F.iso.S),F.Cx(end,end),'x'),legend('isotherme donne par Bi','isotherme fourni','coins de l isotherme fourni','point initial',strcat(['fin simulation a t=',num2str(F.tC(end))])),xlabel('aw'),ylabel('C')
                figure(4)
		clf
		title('verification par dérivée numérique de l equation de Fick : haut et bas doivent coincider.')
		subplot(224),hold on,for not=1:maxnot;plot(ddCxsurdxdx(not,:),F.x(2:end-1),'color',tcolor(not,:)),end,ylabel('x'),xlabel('ddCsurdxdx')
		subplot(223),hold on,for nox=2:maxnox-1;plot(respde.tC,ddCxsurdxdx(:,nox-1),'color',xcolor(nox,:)),end,xlabel(Fostring),ylabel('ddCsurdxdx')
		subplot(222),hold on,for not=2:maxnot-1;plot(dCxsurdt(not-1,:),F.x,'color',tcolor(not,:)),end,ylabel('x'),xlabel('dCsurdt')
		subplot(221),hold on,for nox=1:maxnox;plot(respde.tC(2:end-1),dCxsurdt(:,nox),'color',xcolor(nox,:)),end,xlabel(Fostring),ylabel('dCsurdt')
	end
	if any(strcmp(lower(methode),{'all' 'odepde' 'ode'}))
		figure(5)
		subplot(232), hold on, plot(ti,Codei,'r-')
		subplot(233), hold on, plot(sqrt(ti),Codei,'r-')
		subplot(234), hold on, plot(ti,fodei,'r-')
		subplot(236), hold on, plot(Codei,fodei,'m-'), xlabel('C'), ylabel('j')
		figure(6)
		subplot(121), hold on, plot(ti/F.Bi,Codei,'k-')
		subplot(122), hold on, plot(Codei,1/F.Bi*fodei,'k-')
	end
	if any(strcmp(lower(methode),{'all' 'odepde'}))
		figure(7)
		subplot(235), hold on, plot(Cmi,Codei,'ro',[min(Cmi) 1],[min(Codei) 1],'k-')
		xlabel('C_{pde}'), ylabel('C_{ode}')
	end
	if any(strcmp(lower(methode),{'all' 'analytic'}))
		figure(8)
		clf
		subplot(231), hold on, plot(F.x*ones(1,size(resanalytic.tC,1)),resanalytic.Cx'), xlabel('x'), ylabel('C')
		subplot(232), hold on, plot(resanalytic.tC(1:end-1),resanalytic.f4tC(1:end-1),'b-'), ylabel('j'),xlabel(Fostring)
		subplot(233), hold on, plot(F.C0sol-resanalytic.fc4tC/F.adimensionne.L1,resanalytic.f4tC,'b-'), xlabel('Cmoyen'), ylabel('j'), title('kpd')
		subplot(234), hold on, plot(sqrt(resanalytic.tC(1:end-1)),F.C0sol-resanalytic.fc4tC(1:end-1)/F.adimensionne.L1,'b-'), xlabel([Fostring,'^{1/2}']),ylabel('Cmoyen')
		subplot(235), hold on, plot(resanalytic.tC(1:end-1),F.C0sol-resanalytic.fc4tC(1:end-1)/F.adimensionne.L1,'b-'), xlabel(Fostring),ylabel('Cmoyen')
		subplot(4,3,12), hold on, plot(F.C0sol-resanalytic.fc4tC/F.adimensionne.L1,(-dddsenscamsurdtdtdt(resanalytic.tC,resanalytic.parametres)-ddsenscamsurdtdt(resanalytic.tC,resanalytic.parametres).^2+eps^4)./(eps^4+resanalytic.f4tC.^2),'b-'), ylabel('d2j/dCmoyen2')
		subplot(4,3,9), hold on, plot(F.C0sol-resanalytic.fc4tC/F.adimensionne.L1,(eps^4+ddsenscamsurdtdt(resanalytic.tC,resanalytic.parametres))./(eps^4+resanalytic.f4tC),'b-'), ylabel('dj/dCmoyen')
		figure(9)
		clf
		title('kpd'),hold on
		subplot(121), hold on, plot(resanalytic.tC(1:end-1)/F.Bi,resanalytic.Cm4tC(1:end-1),'k:'), xlabel([Fostring,'/Bi']), ylabel('C')
		subplot(122), hold on, plot(resanalytic.Cm4tC,1/F.Bi*resanalytic.f4tC,'k:'), xlabel('C'), ylabel('j/Bi')
		figure(10)
		clf
		title('batterie de vérifications numériques.'),hold on
		maxnot=size(resanalytic.tC,1),
		maxnox=size(F.x,1),
		tcolor=interp1(0:11,jetcolormap12(0),11*(0:maxnot-1)/(maxnot-.9),'linear');
		xcolor=interp1(0:11,jetcolormap12(0),11*(0:maxnox-1)/(maxnox-.9),'linear');
		subplot(222),hold on,for not=1:maxnot;plot(F.x,resanalytic.Cx(not,:),'color',tcolor(not,:)),end,xlabel('x'),ylabel('C'),title('la dérivée en 0 doit entre nulle')
		subplot(221),hold on,for nox=1:maxnox;plot(resanalytic.tC(1:end-1),resanalytic.Cx((1:end-1),nox),'color',xcolor(nox,:)),end,xlabel(Fostring),ylabel('C'),title('il doit y avoir continuité aux temps de coins d isothermes')
		%i;// ensia ensbana initia agro : veronique a torena ?
		Cxpourdt=[resanalytic.Cx,resanalytic.tC.^2/2];
		dCxsurdt=spdiags((resanalytic.tC(3:end)-resanalytic.tC(1:end-2)).^(-1),0,maxnot-2,maxnot-2)*(Cxpourdt(3:end,:)-Cxpourdt(1:end-2,:));
		%i;// pour que interp1 n'ait pas deux fois la meme abscisse.
		trouvezeros=find(diff(dCxsurdt(:,end))<0);
		dCxsurdt(trouvezeros+1,end)=min(dCxsurdt(trouvezeros,end)+eps^.8*max(1,abs(dCxsurdt(trouvezeros,end))),(dCxsurdt(trouvezeros+2,end)+dCxsurdt(trouvezeros,end))/2);
		try
		  dCxsurdt=interp1(dCxsurdt(:,end),dCxsurdt(:,1:end-1),resanalytic.tC(2:end-1),'linear');
		  Cxpourdxdx=[resanalytic.Cx;F.x'.^3/6];
		  dCxsurdxgauche=(Cxpourdxdx(:,2:end)-Cxpourdxdx(:,1:end-1))*spdiags((F.x(2:end)-F.x(1:end-1)).^(-1),0,maxnox-1,maxnox-1); %i;// dCxsurdx decentre a droite
		  ddCxsurdxdx=(dCxsurdxgauche(:,2:end)-dCxsurdxgauche(:,1:end-1))*spdiags((F.x(3:end)-F.x(1:end-2)).^(-1),0,maxnox-2,maxnox-2)*2; %i;// et on redérive. Ce n'est pas tellement centré :-(
		  ddCxsurdxdx=interp1(ddCxsurdxdx(end,:),ddCxsurdxdx(1:end-1,:)',F.x(2:end-1),'linear')';
		  errrelative=(ddCxsurdxdx(2:end-1,:)*F.adimensionne.Diff-dCxsurdt(:,2:end-1))./(eps^.5+max(abs(ddCxsurdxdx(2:end-1,:))*F.adimensionne.Diff,abs(dCxsurdt(:,2:end-1))));
		  errrelative=max(0,10+log10(eps+abs(errrelative))).*sign(errrelative);
		  %i;// listetrait=kronenlogs([.2 .3 .5 .7],0*[-6:0]).*(10.^kronenlogs(0*[.2 .3 .5 .7],[-6:0]));
		  subplot(223),hold on;contour(resanalytic.tC(2:end-1)',F.x(2:end-1),errrelative',[-9.75:.5:10]);contour(resanalytic.tC(2:end-1)',F.x(2:end-1),errrelative',[-9.5:1:10]);[h,c]=contour(resanalytic.tC(2:end-1)',F.x(2:end-1),errrelative',[-10:10]);clabel(h,c);xlabel(Fostring),ylabel('x'), title('(10+log10)*sign residu relatif : doit etre abs<7(bruit differentiation numerique de Cx) sauf la ou residu absolu <eps^.25')
		  for not=1:maxnot-1;plot(resanalytic.tC(not),max(F.x),'color',tcolor(not,:));end
		  for nox=1:maxnox;plot(max(resanalytic.tC(1:end-1)),F.x(nox),'color',xcolor(nox,:));end
		  errabsolue=(ddCxsurdxdx(2:end-1,:)*F.adimensionne.Diff-dCxsurdt(:,2:end-1));
		  errabsolue=log10(eps+abs(errabsolue));
		  subplot(4,4,16),hold on;contour(resanalytic.tC(2:end-1)',F.x(2:end-1),errabsolue',[-19.75:.5:20]);contour(resanalytic.tC(2:end-1)',F.x(2:end-1),errabsolue',[-19.5:1:20]);[h,c]=contour(resanalytic.tC(2:end-1)',F.x(2:end-1),errabsolue',[-20:20]);clabel(h,c);xlabel(Fostring),ylabel('x'), title('log10(abs(residu absolu)) calculé par Cx')
		  for not=1:maxnot-1;plot(resanalytic.tC(not),max(F.x),'color',tcolor(not,:));end
		  for nox=1:maxnox;plot(max(resanalytic.tC(1:end-1)),F.x(nox),'color',xcolor(nox,:));end
		subplot(4,4,11),fluxbord=(2*F.x(end)-F.x(end-1)-F.x(end-2))'*dCxsurdxgauche(1:end-1,end)+(F.x(end)-F.x(end-1))'*dCxsurdxgauche(1:end-1,end-1);
		fluxbord=-F.adimensionne.Diff*fluxbord./(2*F.x(end)-F.x(end-1)-F.x(end-2)+F.x(end)-F.x(end-1));
		fluxbordbis=-(resanalytic.fcb4tC(3:end)-resanalytic.fcb4tC(1:end-2))./(resanalytic.tC(3:end)-resanalytic.tC(1:end-2));
		fluxbordter=(resanalytic.fc(3:end)-resanalytic.fc(1:end-2))./(resanalytic.t(3:end)-resanalytic.t(1:end-2));
		fluxbordquart=resanalytic.f4tC;
		hold on,plot(resanalytic.tC(1:end-1),fluxbord(1:end-1),resanalytic.tC(2:end-1),fluxbordbis,resanalytic.t(2:end-1),fluxbordter,resanalytic.tC(1:end-1),fluxbordquart(1:end-1),resanalytic.tC(2:end-1),fluxbord(2:end-1)-fluxbordbis,resanalytic.t(2:end-1),interp1(resanalytic.tC,fluxbord,resanalytic.t(2:end-1),'linear')-fluxbordter,resanalytic.t(2:end-1),interp1(resanalytic.tC,fluxbordquart,resanalytic.t(2:end-1),'linear')-fluxbordter),legend('flux par Cx','flux par -fcb4tC','flux par fc','flux par f4tC','difference Cx a -fcb4tC','difference Cx a fc','difference f4tC a fc'),xlabel(Fostring),ylabel('flux'),title('les differences doivent etre nulles a eps^.5 pres')
		danscolonne=1;%i;danscolonne='r';// pour scilab
		subplot(4,4,12),hold on,plot(resanalytic.tC(1:end-1),resanalytic.Cm4tC(1:end-1),resanalytic.t,resanalytic.C,resanalytic.t,F.C0sol-resanalytic.fc/F.adimensionne.L1,resanalytic.t(1:end-1),interp1(resanalytic.tC,resanalytic.Cm4tC,resanalytic.t(1:end-1),'linear')-F.C0sol+resanalytic.fc(1:end-1)/F.adimensionne.L1,resanalytic.t(1:end-1),resanalytic.C(1:end-1)-F.C0sol+resanalytic.fc(1:end-1)/F.adimensionne.L1),legend('Cmoyen par Cm4tC','Cmoyen par C','Cmoyen par fc','difference Cm4tC a fc','difference C a fc'),xlabel(Fostring),ylabel('concentration'),title('les differences doivent etre nulles a eps^.5 pres')
		concentrationinterfaceliquide=F.adimensionne.C0flu-resanalytic.Cm4tC*F.L+F.C0sol*F.L+fluxbordquart/F.Bi;
		Cs=[F.iso.S.Cs,resanalytic.parametres(1).a0p]; if Cs(end)<Cs(end-1); Cs=Cs(1:end-1); end;
		subplot(4,4,15),hold on,plot(F.iso.F.invconc(concentrationinterfaceliquide,F.iso.F),resanalytic.Cx(:,end),F.iso.S.invconc(Cs,F.iso.S),Cs,F.iso.S.invconc(Cs,F.iso.S),Cs,'d',F.iso.S.invconc(Cs,F.iso.S),F.iso.S.conc(F.iso.S.invconc(Cs,F.iso.S),F.iso.S),'d',F.iso.F.invconc(F.adimensionne.C0flu,F.iso.F),F.C0sol,'+',F.iso.S.invconc(resanalytic.Cx(end,end),F.iso.S),resanalytic.Cx(end,end),'x',F.iso.F.invconc(concentrationinterfaceliquide,F.iso.F),F.iso.S.conc(F.iso.F.invconc(concentrationinterfaceliquide(:)',F.iso.F),F.iso.S)-resanalytic.Cx(:,end)'),if figreset==true; legend('isotherme par Bi, C0flu, Cm4tC et f4tC','isotherme fourni','coins de l isotherme avec t','verification isotherme inverse','point initial',strcat(['fin simulation a t=',num2str(resanalytic.tC(end))]),'difference par Bi et fourni'),end,xlabel('aw'),ylabel('C'),title('l isotherme doit etre suivi a eps^.5 pres');
		if length(resanalytic.parametres)>2
		  text(F.iso.S.invconc([resanalytic.parametres(3:end).csp],F.iso.S),[resanalytic.parametres(3:end).csp],num2str([resanalytic.parametres(3:end).tp]'));
		end
                figure(11)
		clf
		title('allure des deux membres de l equation de Fick donne par Cx : haut et bas doivent coincider.'),hold on
		subplot(224),hold on,for not=1:maxnot;plot(ddCxsurdxdx(not,:),F.x(2:end-1),'color',tcolor(not,:)),end,ylabel('x'),xlabel('ddCsurdxdx')
		subplot(223),hold on,for nox=2:maxnox-1;plot(resanalytic.tC(1:end-1),ddCxsurdxdx(1:end-1,nox-1),'color',xcolor(nox,:)),end,xlabel(Fostring),ylabel('ddCsurdxdx')
		subplot(222),hold on,for not=2:maxnot-1;plot(dCxsurdt(not-1,:),F.x,'color',tcolor(not,:)),end,ylabel('x'),xlabel('dCsurdt')
		subplot(221),hold on,for nox=1:maxnox;plot(resanalytic.tC(2:end-1),dCxsurdt(:,nox),'color',xcolor(nox,:)),end,xlabel(Fostring),ylabel('dCsurdt')
		catch
		end
		drawnow
		if dispon
		  %i;// keyboard
		end
	end
end
%i;endfunction
%i;if isdef('%')==%f;function %(varargin);
%i;endfunction;end
%i; for nom=['sensca' 'senscam' 'dsenscamsurdt' 'mesh1D' 'kronenlogs']; monfich=dir(nom+'.m'); toload=%f; if isdef(nom)==%f; toload=%t; else; execstr('if monfich.date~=lastload.'+nom+';toload=%t;end'); end; if toload==%t; execstr('lastload.'+nom+'=monfich.date'); exec(nom+'.m',-1); end; end; // autoload of functions in scilab, like matlab ...

%i;// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           SUB-FUNCTIONS
%i;// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numerical differentiation
% ordre1 [h6y(8)/140]
% ordre2 [h6y(8)/560]
function dydt = ndf(t,y,ordre,dydt0)
%i;[nargout,nargin]=argn(0);// compatibilité scilab.
if nargin<3, ordre = 1; end
if nargin<4, dydt0 = []; end
dt	= t(2)-t(1);
switch ordre;%i;select ordre;// compatibilité scilab matlab.
case 1
	D	= [-1 +9 -45 0 +45 -9 +1]'/(60*dt);
	D0	= [-137 300 -300 200 -75 12]'/(60*dt);
	D1	= -flipud(D0);
case 2
	D	= [2 -27 +270 -490 270 -27 +2]'/(180*dt.^2);
	D0	= [45 -154 +214 -156 61 -10]'/(12*dt.^2);
	D1 = flipud(D0);
end
%i;// yfull	= [repmat(y(1),3,1) ; y ; repmat(y(end),3,1) ];
%i;// dydt	= [yfull(1:end-6) yfull(2:end-5) yfull(3:end-4) yfull(4:end-3) yfull(5:end-2) yfull(6:end-1)  yfull(7:end)] * D;
dydt = 	[
		[y(1:3) y(2:4) y(3:5) y(4:6) y(5:7) y(6:8) ] * D0;
		[y(1:end-6) y(2:end-5) y(3:end-4) y(4:end-3) y(5:end-2) y(6:end-1)  y(7:end)] * D
		[y(end-7:end-5) y(end-6:end-4) y(end-5:end-3) y(end-4:end-2) y(end-3:end-1)  y(end-2:end)] * D1
	];
if any(dydt0), dydt(1) = dydt0; end
dy = diff(dydt);
i = intersect(find(dy*sign(mean(sign(dy)))<0),2:10);
if any(i), ni = setdiff(1:length(t),i); dydt(i) = interp1(t(ni),dydt(ni),t(i),'cubic','extrap'); end
%i;//  dydt 	= interp1(t(6:end-5),dydt(6:end-5),t,'cubic','extrap');
% if mean(diff(dydt(end-3:end)))*mean(diff(dydt(5:end-4)))<0
%i;// 	disp('**')
% end
%i;endfunction

% ODE
function dCdt = odediff(t,C,S)
global timeout_ode
isvectorized = 0;
if toc>timeout_ode, evalc('erron'), end
if nargin<3
	isvectorized = 1;
	S = C;
	C = t;
end
if ~isvectorized
error 'plus faisable comme le champ K a disparu'
	j1		= S.Bi * ( (S.K+S.L)*C(1) - S.L*S.C0sol ) / ( 1 + 2/6 * S.Bi * S.K);
	b		= S.Bi * S.K * sqrt( 3/2 * (S.C0sol - C(1) ) );
	c		= S.Bi * ( (S.K-S.L) * S.C0sol + S.L * C(1) );
	d		= b^2 + 4 * c;
    j2		= (( -b + sqrt(d) ) / 2)^2;
    if j2<=0
        dCdt = -j1;
    else
        if S.C0sol-sqrt(3/2*(S.C0sol-C)/j2)>S.xcrit,		dCdt = -j2;
        else,										dCdt = -j1; end
    end
else
	j1		= S.Bi * ( (S.K+S.L)*C - S.L*S.C0sol ) / ( 1 + 2/6 * S.Bi * S.K);
	b		= S.Bi * S.K * sqrt( 3/2 * (S.C0sol - C ) );
	c		= S.Bi * ( (S.K-S.L) * S.C0sol + S.L * C );
	d		= b.^2 + 4 * c;
	j2		= (( -b + sqrt(d) ) / 2).^2;
    j2(j2<=0) = NaN;
	i2		= S.C0sol-sqrt(3/2*(S.C0sol-C)./j2)>S.xcrit;
	i1		= S.C0sol-sqrt(3/2*(S.C0sol-C)./j2)<=S.xcrit;
	dCdt	= zeros(size(C));
	dCdt(i2) = -j2(i2);
	dCdt(i1) = -j1(i1);
end
%i;endfunction


%i;// PDE problem :  [C,F,S] = PDEFUN(X,T,U,DUDX,P1,P2...)
function [c,f,s] = pdediff_PDE(x,t,U,dUdx,S)
global timeout_pde
if toc>timeout_pde, evalc('erron'), end
c = 1;
f = dUdx;
s = 0;
%i;endfunction

%i;// IC : U = ICFUN(X,P1,P2...)
function U0 = pdediff_IC(x,S)
% concentration initiale dans le solide.
U0 = ones(size(x))*S.C0sol;
%i;endfunction

%i;// BC : [PL,QL,PR,QR] = BCFUN(XL,UL,XR,UR,T,P1,P2...)
function [pa,qa,pb,qb,purge] = pdediff_BC(xa,Ua,xb,Ub,t,Fc,S)
global Hsmooth
% pervious contact on the left
% pb = 0;
% qb = 1; ordre de la condition aux limites.
% pa = S.Bi*( S.K*(Ua-0) - S.L * -Fc(1) );
% qa = -1; ordre de la condition aux limites.
% impervious contact on the left
pa = 0;
qa = 1;
% calc aw on the right
awb=S.iso.S.invconc(Ub,S.iso.S);
CLi = S.iso.F.conc(awb,S.iso.F.PvsatsurRT); % concentration dans le fluide a la frontiere
pb = S.Bi*( (CLi-0) - S.L * -Fc(2) );   
qb = 1;
%i;endfunction

%i;// General 1D PDE formulation (without advection)
%i;//
%i;//  c(x,t,u,Du/Dx) * Du/Dt = x^(-m) * D(x^m * f(x,t,u,Du/Dx))/Dx + s(x,t,u,Du/Dx)
%i;//  --------------                            --------------       --------------
%i;//   capacitance                                   flux                source
%i;// (diagonal positive nxn matrix)                nx1 vector           nx1 vector
%i;//
%i;//  u(t,x) can be defined in R^n on a<=x<=b during t0 <= t <= tf
%i;//  m = 0, 1, or 2 => {slab, cylindrical, or spherical} symmetry
%i;//  subjected to BC:  p(x,t,u) + q(x,t) * f(x,t,u,Du/Dx) = 0
%i;//                               -------
%i;//                           diagonal nxn matrix
%i;//
%i;//  vim:softtabstop=2:shiftwidth=2:guioptions+=a;
function debug_plot(varargin);
  %i;// dbstack
  for i=1:length(varargin)
    instr='arg=varargin{i};';%i;instr=strrep(strrep(instr,'{','('),'}',')');// compatibilite scilab
    evalc(instr);
    disp(size(arg))
  end
  %i;// this line is scilab only
  instr='plot(varargin{:})';%i;instr=strrep(strrep(instr,'{','('),'}',')');// compatibilite scilab
  evalc(instr);
%i;endfunction

function [valeur,stop,monotonie]=centiC0(t,y,varargin)
  %to have a good-looking kpd graph ...
  F=varargin{end};
  % DG 310507 : il faudra comparer ces deux fonctions et voir laquelle est la plus rapide.
  if 1==1;%i;// vrai si on utilise une seule fonction, faux sinon.
    listeC=F.C0sol;
    valeur=cos((y(1,end-2)-listeC)/(F.C0sol-F.iso.S.C_a)*100*pi);
    stop=0;
    monotonie=0;
  else
    listeC=[0.005:0.01:2]'*(F.C0sol-F.iso.S.C_a)+F.iso.S.C_a;
    valeur=y(1,end-2)-listeC;
    stop=0*listeC;
    monotonie=1+0*listeC;
  end
%i;endfunction
function res=xplusltanxmoinsbxxtanx(x,LsurKp,unsurBiKp,n)
  res=x+(LsurKp-unsurBiKp*x^2)*tan(x);
%i;endfunction
function [res,grad,cossin]=xplusnpicosxplusLsurKpsinxmoinsunsurBIKpcarrexplusnpisinx(x,LsurKp,unsurBiKp,n)
  %i;// this function can also take line vectors.
  res=[x+n*pi;LsurKp-unsurBiKp*(x+n*pi).^2];
  res=res./kron([1;1],(x+n*pi));
  normres=sqrt((1+n*pi).^2+(LsurKp-unsurBiKp*(1+n*pi)^2).^2./(1+n)^2);
  cossin=[cos(x);sin(x)];
  danscolonne=1;%i;danscolonne='r';// pour scilab
  grad=sum(res.*[-cossin(2,:);cossin(1,:)]+[zeros(1,size(x,2));-LsurKp./(x+n*pi).^2-unsurBiKp].*cossin,danscolonne)./normres;
  res=sum(res.*cossin,danscolonne)./normres;
%i;endfunction
function res=A_n0__p(lambda_np,lambda_np_cache,Bi,Kp,L);
  res=L./(lambda_np.*lambda_np_cache.sin/Bi+Kp.*lambda_np.*(.5./lambda_np_cache.sin-lambda_np_cache.cos./lambda_np./2));
  %i;// NB: c'est normal qu'on ne puisse pas passer à la limite pour lambda qui tend vers 0.
%i;endfunction
function res=A_n0__psurL(lambda_np,lambda_np_cache,Bi,Kp,L);
  res=1./(lambda_np.*lambda_np_cache.sin/Bi+Kp.*lambda_np.*(.5./lambda_np_cache.sin-lambda_np_cache.cos./lambda_np./2));
  %i;// NB: c'est normal qu'on ne puisse pas passer à la limite pour lambda qui tend vers 0.
%i;endfunction
function [res,lambdan0_cache,minpointfixe]=lambda_n__p(Bi,Kp,L,nmax,varargin)
%i;// calcule la liste des valeurs propres.
  ploton=false;
  first=1;%i;first=[]
  if length(varargin)>0 ;arg=varargin(1); ploton=[arg{first}]; end;
  if ploton;figure(15); clf;end;
  lambdan0MODpi=ones(nmax+2,1);
  lambdan0DIVpi=[0:nmax+1]';
  %i;// lambdan0moinsdemipiDIVpi=[0:nmax+1]';
  for i=1:10;
    lastlambdan0MODpi=lambdan0MODpi;
    %i;// avec acotan, ca aurait été plus stable, notamment car les lambdas sont tous entre 0 et pi, mais TODO comme alors lambdan0MODpi tend vers0, ca aurait été toucher les bornes maximums de cotan, et alors que peut se passer ??? ah, si j'avais une fonction atancotan, ce serait l'idéal. Ah j'ai trouvé, c'est atan2(r*sin(t),r*cos(t)) qui renvoie t pour r>0 et -pi<t<pi.
    lambdan0MODpi=atan2(Bi*Kp*(lambdan0DIVpi*pi+lambdan0MODpi),(lambdan0DIVpi*pi+lambdan0MODpi).^2-Bi*L);
    if any(lambdan0MODpi<0)
      error('erreur bug35 dans lambda_n__p, lambdan0MODpi aurait du rester positif car y est positif dans atan2(y,x)')
      lambdan0MODpi(find(lambdan0MODpi<0))=lambdan0MODpi+pi;
    end
    %i;// lambdan0MODpi=atan(Bi*Kp*(lambdan0DIVpi*pi+lambdan0MODpi)./((lambdan0DIVpi*pi+lambdan0MODpi).^2-Bi*L));
    %i;// lambdan0MODpi(find(lambdan0MODpi<0))=lambdan0MODpi+pi;
    if ploton
      subplot(311)
      plot (log10(1:nmax+2),log10(eps^10+(eps^5+abs(lambdan0MODpi-lastlambdan0MODpi))./abs(lambdan0MODpi)))
      hold on
      subplot(312)
      plot (log10(1:nmax+2),lambdan0MODpi)
      hold on
      subplot(313)
      plot (log10(1:nmax+1),diff(lambdan0MODpi))
      hold on
      drawnow
    end
  end
  if ploton
    subplot(311)
  end
  if Kp<0 | L<0 | Bi<0
    error 'bug0dans lambda_n_'
  end
  %i;// matlab7 optimise et transforme minpointfixe=min(nmax+2,max(find(abs(lastlambdan0MODpi-lambdan0MODpi)>eps^.9*abs(lambdan0MODpi)))+5) en minpointfixe=min(nmax+2,max(find(abs(lastlambdan0MODpi-lambdan0MODpi)>eps^.9*abs(lambdan0MODpi))+[])) de temps en temps ...
  minpointfixe=max(find(abs(lastlambdan0MODpi-lambdan0MODpi)>eps^.9*abs(lambdan0MODpi)));
  %i;// lambdan0MODpi(find(lambdan0MODpi<0))=lambdan0MODpi(find(lambdan0MODpi<0))+pi;
  if isempty(minpointfixe)
    minpointfixe=0;
  end
  minpointfixe=min(nmax+2,minpointfixe+5);
  arefaire=1:minpointfixe;
%i;//   lambdaplusdemipiDIVpiarefaire(2).s=[lambdaplusdemipiDIVpiarefaire(1).s(find((lambdan0DIVpi(lambdaplusdemipiDIVpiarefaire(1).s)*pi+pi/2+3*pi).^2>Bi*L));unique(max(1,[0:3]+floor(sqrt(Bi/L)/pi))];
%i;//   lambdaplusdemipiDIVpiarefaire(1).s=[lambdaplusdemipiDIVpiarefaire(1).s(find((lambdan0DIVpi(lambdaplusdemipiDIVpiarefaire(1).s)*pi-pi/2-3*pi).^2<Bi*L));unique(max(1,[-1:2]+floor(sqrt(Bi/L)/pi))];
  for nolambdan0=arefaire;
    minerr=eps.^(-5);
    optres=-1;
    maxapresdemipi=0+(lambdan0DIVpi(nolambdan0)<ceil(sqrt(Bi*L)/pi));
    apresdemipi=min(nolambdan0*4-7,1-(lambdan0DIVpi(nolambdan0)+1>floor(sqrt(Bi*L)/pi)));
    basex=[10.^[log10(eps):.1:log10(pi)],pi];
%i;//     liste0=0;
%i;//     xtotlist=0;
%i;//     xlist=pi;
%i;//     if nolambdan0>1
%i;//       xlist=[0 pi];
%i;//     end
%i;//     ytotlist=[];
%i;//     maxit=10;
%i;//     while minerr>eps^.8
%i;// 	xlist=sort([xlist,kronenlogs(liste0,[-basex(end:-1:1),0,basex])]);
%i;// 	xlist(find(xlist<0|xlist>pi))=0;
%i;// 	n=lambdan0DIVpi(nolambdan0)+0*xlist;
%i;// 	[ylist,dysurdx,cossin]=xplusnpicosxplusLsurKpsinxmoinsunsurBIKpcarrexplusnpisinx(xlist,L/Kp,1/Bi/Kp,lambdan0DIVpi(nolambdan0));
%i;// 	liste0=contourc(log10(xlist),[1 2],[1 1]'*ylist,[0 0]);
%i;// 	liste0=[lambdan0MODpi(1),lastlambdan0MODpi(1),10.^liste0(1,2:3:end)];
    dodisplay='off';if ploton;dodisplay='on';end;
    %i;// TODO, attention avec matlab v6, fzero est tres tres lent ...
    [optres,minerr]=fzero(@xplusnpicosxplusLsurKpsinxmoinsunsurBIKpcarrexplusnpisinx,[eps,pi],optimset('Display',dodisplay,'TolX',eps^5),L/Kp,1/Bi/Kp,lambdan0DIVpi(nolambdan0));
    apresdemipi=maxapresdemipi+2;
    while apresdemipi<=maxapresdemipi
      if apresdemipi==-3
	%i;// [1,x].*[x+n*pi;LsurKp-unsurBiKp*(x+n*pi).^2];
	xlist=[10.^[log10(eps):.1:log10(pi)],pi];
	n=lambdan0DIVpi(nolambdan0)+0*xlist;
	[ylist,dysurdx,cossin]=xplusnpicosxplusLsurKpsinxmoinsunsurBIKpcarrexplusnpisinx(xlist,L/Kp,1/Bi/Kp,lambdan0DIVpi(nolambdan0));
	liste0=contourc(log10(xlist),[1 2],[1 1]'*ylist,[0 0]);
	liste0=[lambdan0MODpi(1),lastlambdan0MODpi(1),10.^liste0(1,2:3:end)];
	if ploton
	  figure(24);
	  subplot(221)
	  plot(xlist,(L/Kp-1/Bi/Kp*(xlist+n*pi).^2)./(xlist+n*pi));
	  subplot(222)
	  plot(xlist,dysurdx)
	  subplot(223)
	  plot(xlist,cossin(1,:)./cossin(2,:));
	  subplot(224)
	  hold on
	  plot(xlist,ylist)
	  plot(liste0,0*liste0,'+');
	end
	etendu=0;
	%i;// keyboard
	apresdemipi=-2;
      end
      if apresdemipi==-2;
	if isempty(liste0)
	  apresdemipi=-1;
	else
	  xinitial=liste0(1);
	  if etendu==1
	    lambdamin=0;
	    lambdamax=pi/2;
	  else
	    lambdamin=max([0,xlist(find(xlist<xinitial))]);
	    lambdamax=min([pi,xlist(find(xlist>xinitial))]);
	    liste0=liste0(2:end);
	  end
	end
      end
      if apresdemipi>-2;
	lambdamin=max(eps,apresdemipi*pi/2);
	lambdamax=max(1,(1+apresdemipi))*pi/2;
	xinitial=lambdan0MODpi(nolambdan0);
	apresdemipi=apresdemipi+1;
      end
      if apresdemipi==0;
	xinitial=sqrt((L/Kp-1)*Bi*Kp)-lambdan0DIVpi(nolambdan0)*pi;
	%i;// developpement limite de [1,x].*[1;LsurKp-unsurBiKp*(x+n*pi).^2];
      end
  %i;//   lambdaplusdemipiDIVpiarefaire(1).s=[lambdaplusdemipiDIVpiarefaire(1).s(find((lambdan0DIVpi(lambdaplusdemipiDIVpiarefaire(1).s)*pi-pi/2-3*pi).^2<Bi*L));unique(max(1,[-1:2]+floor(sqrt(Bi/L)/pi))];
      %i;// L du livre = alphadiametresurD/2
      err=1;
      maxit=1;
      while (err>eps^.8) & (maxit>0)
%i;//       lamin=lambdamin;
%i;//       %i;// specifique pour ne pas qu'il trouve la solution double x=0 ... et ca ne sert pas toujours
%i;//       if i==1
%i;// 	if L/Kp-1/Bi/Kp*(pi+(i-1)*pi)^2>0
%i;//           lamin=pi/2;
%i;// 	end
%i;//       end
	%i;// try
	  [maxXmh,err]=fzero(@xplusnpicosxplusLsurKpsinxmoinsunsurBIKpcarrexplusnpisinx,[lambdamin,lambdamax],optimset('Display',dodisplay,'TolX',eps^5),L/Kp,1/Bi/Kp,lambdan0DIVpi(nolambdan0))
%i;//	  [maxXmh,err]=lsqcurvefit(@xplusnpicosxplusLsurKpsinxmoinsunsurBIKpcarrexplusnpisinx,xinitial,L/Kp,0,lambdamin,lambdamax,optimset('Display',dodisplay,'Jacobian','on','DerivativeCheck','off'),1/Bi/Kp,lambdan0DIVpi(nolambdan0));
%i;//	  if ploton
%i;//	    figure(24);
%i;//	    subplot(221)
%i;//	    plot([xinitial,maxXmh],[apresdemipi,err],'r')
%i;//	    disp([apresdemipi,1,xinitial,maxXmh,err,log10(err)]);
%i;//	  end
%i;//	  [maxXmh,err]=lsqcurvefit(@xplusnpicosxplusLsurKpsinxmoinsunsurBIKpcarrexplusnpisinx,maxXmh,L/Kp,0,lambdamin,lambdamax,optimset('Display',dodisplay,'TolFun',eps^5,'TolX',eps^5,'MaxIter',44,'Jacobian','on','DerivativeCheck',dodisplay),1/Bi/Kp,lambdan0DIVpi(nolambdan0));
	  if ploton
	    figure(24);
	    subplot(221)
	    plot([xinitial,maxXmh],[apresdemipi,err],'g')
	    disp([apresdemipi,2,xinitial,maxXmh,err,log10(err)]);
	  end
	  if minerr>abs(err)
	    minerr=abs(err);
	    if ploton
	      figure(24);
	      subplot(221)
	      plot([optres,maxXmh],[apresdemipi,err],'k')
	    end
	    optres=maxXmh;
	  end
	%i;// catch
	%i;// end
	%i;// [maxXmh,err]=fsolve(@xplusltanxmoinsbxxtanx,xinitial,optimset('Display',DoDisplay),[L/Kp,1/Bi/Kp]); %i;// fsolve ne convient pas car il me faut des bornes.
	maxit=maxit-1;
	xinitial=rand(1,1)*(lambdamax-lambdamin)+lambdamin;
      end
      %i;// TODO : faire une version rapide à trois points, on calcule un point intermediare soit le milieu soit sécante, puis on fait un pas de gradient conjugué, puis on recommence, cette version étant rapide car elle fait tout vectoriellement, pour tous les nolambdan0 à la fois.
      %i;// TODO ou faire encore mieux, avec lambdan0MODpi=arctan(fractionde{\[ \tan (\lambda ) = \frac{{{\text{Bi}}K_1 \lambda }} {{\lambda ^2  - {\text{Bi}}L}} \]}([0:nmax+1]+lambdan0MODpi)) repete deux ou trois fois, en partant de 0, et le tour est joue !
    end
    if optres==-1
      error 'bugdans lambda_n_'
      maxXmh=-1;
    end
    if ploton
      maxXmh=optres;
      figure(15)
      plot(log10(nolambdan0),log10(eps^4+abs(lambdan0MODpi(nolambdan0)-maxXmh)),'+');
      plot(log10(nolambdan0),log10(eps^3+abs(maxXmh)),'x');
    end
    lambdan0MODpi(nolambdan0)=optres;
%i;//    if maxXmh==lambdamin
%i;//      error 'bug4dans solution double pour nolambdan0==1 dans lambda_n_'
%i;//    end
  end
  lambdan0=lambdan0MODpi+[0:nmax+1]'*pi;
  res=unique(lambdan0(find(lambdan0>0)));
  if max(abs(size(res)-size(lambdan0)))>0
    error 'bug11dans lambda_n_ size'
  end
  if max(abs(res-lambdan0))>0
    error 'bug12dans lambda_n_ sort'
  end
  res=res(1:nmax);
  lambdan0_cache.modpi=lambdan0MODpi(1:nmax);
  lambdan0_cache.mod2pi=lambdan0_cache.modpi+mod([0:nmax-1]',2)*pi;
  lambdan0_cache.sin=sin(lambdan0_cache.mod2pi);
  lambdan0_cache.cos=cos(lambdan0_cache.mod2pi);
  lambdan0_cache.energy=res.^2.*lambdan0_cache.sin.^2+Bi*Kp*res.^2.*(1-(eps^5+lambdan0_cache.sin)./(eps^5+res).*lambdan0_cache.cos)/2;
  if max(diff(res))>1.5*pi
    error 'bug2dans lambda_n_'
  end
  if max(abs(size(lambdan0_cache.modpi)-size(res))>0)
    error 'bug3dans lambda_n_'
  end
%i;endfunction
function res=inlinecmp(a,b)
  res=strcmp(char(a),char(b))&all(strcmp(argnames(a),argnames(b)));
%i;endfunction
function [res,gradient]=trouvecoin(logt,segment,cspplus1)
  res=sensca(exp(logt)+segment.tp,segment.x,segment)-cspplus1;
  %i;// gradient=gradienttbordanalytic(exp(logt),segment).*exp(logt)
%i;endfunction
function [res,last20]=gradienttbordanalytic(t,segment)
  last20=cos(segment.lambda_np*segment.x).*(exp(-(t-segment.tp)*segment.lambda_np.^2).*segment.lambda_np.^2);
  res=segment.a_np'*last20;
  %i;// TODO mettre expm1 qui est l'analogue de log1p, mais matlab65 semble le confondre avec l'exponentielle de matrice :-( donc pas facile à utiliser.
%i;endfunction
function [res,last20]=primitivetgradientxbordanalytic(t,segments)
  last20=0;
  t=t(:);
  res=0*t;
  for not=1:size(t,1)
    for p=1:length([segments.tp])
      segment=segments(p);
      if segment.tp<t(not)
	tmax=t(not);
	if p<length([segments.tp])
	  tmax=min(tmax,segments(p+1).tp);
	end
	tmp=segment.lambda_np.*segment.lambda_np_cache.sin.*(1-exp(-(tmax-segment.tp)*segment.lambda_np.^2))./segment.lambda_np.^2;
	res(not)=res(not)+segment.a_np'*tmp;
	%i;// TODO adapter last20 pour quand t a plusieurs coefficients.
	last20=last20+segment.a_np(end-segment.lastterms+1:end).*tmp(end-segment.lastterms+1:end);
      end
    end
  end
  %i;// TODO mettre expm1 qui est l'analogue de log1p, mais matlab65 semble le confondre avec l'exponentielle de matrice :-( donc pas facile à utiliser.
%i;endfunction
function [Cres,Clast20]=sensca(t,x,segments)
%i;// analytic solution evaluated by user.
  t=t(:);
  x=x(:);
  Cres=zeros(size(t,1),size(x,1));
  Clast20=zeros(size(t,1),size(x,1)*segments(1).lastterms);
  for not=1:size(t,1)
    p=max(find(vec2list([segments.tp])<t(not)));
    segment=segments(p);
    for nox=1:size(x,1)
      last20=cos(segment.lambda_np*x(nox)).*exp(-(t(not)-segment.tp)*segment.lambda_np.^2);
      Cres(not,nox)=segment.a0p+segment.a_np'*last20;
      if segment.lastterms>0
        Clast20(not,(nox-1)*segment.lastterms+(1:segment.lastterms))=(segment.a_np(end-segment.lastterms+1:end).*last20(end-segment.lastterms+1:end))';
      end
    end
  end
%i;endfunction
%i;if 1==0
function res=vec2list(vec)
  res=vec;
%i;endfunction;end
function [Cres,Clast20]=senscam(t,segments)
  % SENSC Analytic Moyenné au temps t si segments=resanalytic.parametres renvoye par SENSC('analytic')
  t=t(:);
  Cres=zeros(size(t,1),1);
  Clast20=zeros(size(t,1),1);
  for not=1:size(t,1)
    p=max(find(vec2list([segments.tp])<t(not)));
    segment=segments(p);
    x=segment.x(:);
    for nox=1:size(x,1)
      last20=sin(segment.lambda_np*segment.x(nox))./segment.lambda_np.*exp(-(t(not)-segment.tp)*segment.lambda_np.^2);
      Cres(not,nox)=segment.a0p*x(nox)+segment.a_np'*last20;
      if segment.lastterms>0
	Clast20(not,(nox-1)*segment.lastterms+(1:segment.lastterms))=(segment.a_np(end-segment.lastterms+1:end).*last20(end-segment.lastterms+1:end))';
      end
    end
  end
%i;endfunction
function [Cres,Clast20]=dsenscamsurdt(t,segments)
  % gradient de SENSC Analytic Moyenné au temps t avec segment renvoye par SENSC('analytic')
  t=t(:);
  Cres=zeros(size(t,1),1);
  Clast20=zeros(size(t,1),1);
  for not=1:size(t,1)
    p=max(find(vec2list([segments.tp])<t(not)));
    segment=segments(p);
    x=segment.x(:);
    for nox=1:size(x,1)
      last20=-segment.lambda_np.*sin(segment.lambda_np*x(nox)).*exp(-(t(not)-segment.tp)*segment.lambda_np.^2);
      Cres(not,nox)=segment.a_np'*last20;
      if segment.lastterms>0
        Clast20(not,(nox-1)*segment.lastterms+(1:segment.lastterms))=(segment.a_np(end-segment.lastterms+1:end).*last20(end-segment.lastterms+1:end))';
      end
    end
  end
%i;endfunction
function [Cres,Clast20]=ddsenscamsurdtdt(t,segments)
  % gradient second de SENSC Analytic Moyenné au temps t avec segment renvoye par SENSC('analytic')
  t=t(:);
  Cres=zeros(size(t,1),1);
  Clast20=zeros(size(t,1),1);
  for not=1:size(t,1)
    p=max(find(vec2list([segments.tp])<t(not)));
    segment=segments(p);
    x=segment.x(:);
    for nox=1:size(x,1)
      last20=segment.lambda_np.^3.*sin(segment.lambda_np*x(nox)).*exp(-(t(not)-segment.tp)*segment.lambda_np.^2);
      Cres(not,nox)=segment.a_np'*last20;
      if segment.lastterms>0
        Clast20(not,(nox-1)*segment.lastterms+(1:segment.lastterms))=(segment.a_np(end-segment.lastterms+1:end).*last20(end-segment.lastterms+1:end))';
      end
    end
  end
%i;endfunction
function [An,normu]=analytic_inverse_Asagiv(segment,L,Bi,Bn)
  Kp=segment.Kp;
  %i;// maple('sensc1:=subs(cos(lambda[n*p])=`segment.lambda_np_cache.cos`, lambda[n*p]=`segment.lambda_np` , K[p]=Kp,EQ65);') ; puis remplacement de ( )^(1/2) par sqrt( ), puis on enleve les ` puis on met des .* et des .^ et des ./ puis on met les réels à la fin.
  sensc1 = 1/2*sqrt(2+4*segment.lambda_np_cache.cos.^2./(segment.lambda_np.^2-Bi*L).*(1/2+Bi*L./(segment.lambda_np.^2-Bi*L))*Bi*Kp);
  %i;// maple('sensc2:=subs(cos(lambda[n*p])=`segment.lambda_np_cache.cos`, lambda[n*p]=`segment.lambda_np` , K[p]=Kp,EQ67);')
  sensc2 = 2*Bi*segment.lambda_np_cache.cos./(segment.lambda_np.^2-Bi*L)./sqrt(2+4*segment.lambda_np_cache.cos.^2./(segment.lambda_np.^2-Bi*L).*(1/2+Bi*L./(segment.lambda_np.^2-Bi*L))*Bi*Kp)*sqrt(Kp*L);
  Bn=(Bn./sensc1);
  Bn=Bn+sensc2*(sensc2'*Bn)/(1-sensc2'*sensc2);
  An=(Bn./sensc1);
  normu=min(abs(sensc1))
  normu=min(find(abs(sensc1)==normu));
  normu=[sqrt(sensc2'*sensc2),sensc1(normu),normu];
%i;endfunction
function [Cres,Clast20]=dddsenscamsurdtdtdt(t,segments)
  % gradient second de SENSC Analytic Moyenné au temps t avec segment renvoye par SENSC('analytic')
  t=t(:);
  Cres=zeros(size(t,1),1);
  Clast20=zeros(size(t,1),1);
  for not=1:size(t,1)
    p=max(find(vec2list([segments.tp])<t(not)));
    segment=segments(p);
    x=segment.x(:);
    for nox=1:size(x,1)
      last20=-segment.lambda_np.^5.*sin(segment.lambda_np*x(nox)).*exp(-(t(not)-segment.tp)*segment.lambda_np.^2);
      Cres(not,nox)=segment.a_np'*last20;
      if segment.lastterms>0
        Clast20(not,(nox-1)*segment.lastterms+(1:segment.lastterms))=(segment.a_np(end-segment.lastterms+1:end).*last20(end-segment.lastterms+1:end))';
      end
    end
  end
%i;endfunction
function iso=baldev02isotherms(name,isoinput)
  first=1;%i;first=[]
  iso=isoinput;
  iso.name=name;
  switch lower(name);%i;select lower(name);// compatibilité scilab matlab.
    case 'bet',
      iso.executive='monolayer, for dehydrated food.';
      iso.source_paramconc='Baldev02(2)';
      iso.paramconc=inline('1./((1-aw)*M)-(1/iso.M_m+1./(iso.C*iso.M_m).*(1-aw)./aw)','M','aw','iso');
      iso.conc=inline('1./(1/iso.M_m+1./(iso.C*iso.M_m).*(1-aw)./aw)./(1-aw)','aw','iso');
    case 'smith'
      iso.executive='biopolymers';
      iso.source_conc='Baldev02(3)';
      iso.conc=inline('iso.M_b-iso.M_a*log(1-aw)','aw','iso');
      iso.invconc=inline('1-exp((iso.M_b-M)/iso.M_a)','M','iso');
    case 'halsey'
      iso.executive='multilayer at distance';
      iso.source_conc='Baldev02(5)';
      iso.conc=inline('exp(iso.a+iso.b*log(-log(aw)))','aw','iso');
      iso.invconc=inline('exp(-exp((log(M)-iso.a)/iso.b))','M','iso');
    case 'caurie'
      iso.executive='used once for gelatin';
      iso.source_conc='Baldev02(6)';
      iso.conc=inline('exp(log(iso.A)-iso.r*aw)','aw','iso');
      iso.invconc=inline('(log(iso.A)-log(M))/iso.r','M','iso');
    case 'bradley'
      iso.executive='first layer, once activated, polarize subsequent layers';
      iso.source_conc='Baldev02(8)';
      iso.conc=inline('log(log(1./aw)/iso.K_2)/log(iso.K_1)','aw','iso');
      iso.invconc=inline('1./exp(exp(log(iso.K_2)+M*log(iso.K_1)))','M','iso');
    case 'oswin'
      iso.executive='sigmoid adapted to sorption';
      iso.source_conc='Baldev02(10)';
      iso.conc=inline('exp(log(iso.a)+iso.n*log(eps^5+aw./(1-aw)))','aw','iso');
      iso.derconc=inline('exp(log(iso.a)+iso.n*log(eps^5+aw./(1-aw))).*(1./(eps^5+aw)+1./(1-aw))*iso.n','aw','iso');
      iso.invconc=inline('1-1./(exp((log(eps^5+M)-log(iso.a))/iso.n)+1)','M','iso');
    case 'langmuir'
      iso.executive='langmuir';
      iso.source_conc='Langmuir; J. Am. Chem. Soc. 38, 2221-95 1916';
      iso.conc=inline('iso.gammamax*iso.K*aw./(1+iso.K*aw)','aw','iso');
      iso.derconc=inline('iso.gammamax*(iso.K./(1+iso.K*aw)-iso.K*iso.K.*aw./(1+iso.K*aw).^2)','aw','iso');
      %i;// TODO voir si la formule suivante est juste : iso.derconc=inline('iso.gammamax*(iso.K/(1+iso.K*aw).^2)','aw','iso');
      iso.invconc=inline('(1./(1-M/iso.gammamax)-1)./iso.K','M','iso');
    otherwise;%i;else
      error(strcat([name,' is an unknown isotherm name.']))
  end
  if ~isfield(iso,'invconc')
    %i;// iso.funfsolve=inline('iso.conc(aw,iso)-M','aw','M','iso');
    %i;// iso.invconc=inline('fsolve(iso.funfsolve,mean(iso.awminawmax),[],M,iso)','M','iso');
    iso.funlsqcurvefit=inline('iso.conc(aw,iso)-M','aw','M','iso');
    iso.invconc=inline('lsqcurvefit(iso.funlsqcurvefit,mean(iso.awminawmax(1:2))+0*M,M,0*M,iso.awminawmax(1)+0*M,iso.awminawmax(2)+0*M,optimset(''TolFun'',eps,''TolX'',eps,''Display'',''off''),iso)','M','iso');
  end
%i; endfunction
function isoc=approximate_isotherm(iso,varargin)
  first=1;%i;first=[]
  nbinter=0; if length(varargin)>0;arg=varargin(1);arg=[arg{first}];if ~isempty(arg); nbinter=arg;end;end;
  ploton=false; if length(varargin)>1;arg=varargin(2);arg=[arg{first}];if ~isempty(arg); ploton=arg;end;end;
  if abs(nbinter(1))>0;
    betiso=iso;
    betiso.awminawmax(3)=-1+sign(1e-77+nbinter(2))-1/abs(nbinter(1));
    %i;// si nbinter(2) est différent de 1 ou -1, on translate les segments en suivant la derivée seconde numérique.
    isoc=approximate_isotherm(betiso);
    isoc.approximate_isotherm_nbinter=nbinter(:)';
    aw=betiso.invconc(isoc.Cs,betiso);
    if nbinter(1)<-1
      %i;// on refait tout, translation de la fonction inverse, 
      cs=[betiso.conc(betiso.aw0flu,betiso),betiso.Ceqsol,betiso.C0sol,betiso.conc(betiso.awminawmax(2),betiso)];
      nbs=max(2,ceil(abs(nbinter(1))/diff(cs(2:3))*diff(cs)));
      isoc.Cs=[]
      for nonbs=1:length(nbs)
	isoc.Cs=real([isoc.Cs(1:length(isoc.Cs)-1),cs(nonbs):diff(cs(nonbs:nonbs+1))/nbs(nonbs):cs(nonbs+1)]);
      end
      aw=betiso.invconc(isoc.Cs,betiso);
      if nbinter(2)>0
	%i;// avec en bonus, tangeantes par les points betiso.aw0flu, betiso.C0sol et betiso.Ceqsol.
	awsp=betiso.invconc(isoc.Cs+eps^.3,betiso);
	awsm=betiso.invconc(isoc.Cs-eps^.3,betiso);
	derseconde=(awsp+awsm-2*aw)/eps^.6;
	isoc.aws=real(aw+(abs(nbinter(2))-1)*derseconde.*([0,diff(isoc.Cs)]+[diff(isoc.Cs),0]).^2);
	disp(isoc.aws);
	%i;// a droite des points de référence :
	ind=cumsum([1,nbs(1:2)]);
	derprim=1./real(betiso.derconc(isoc.aws(ind),betiso));
	df=[diff([isoc.Cs(ind);isoc.Cs(ind+1);isoc.Cs(ind+2)]);diff([isoc.aws(ind);isoc.aws(ind+1);isoc.aws(ind+2)])];
	dCs=real((df(1,:).*df(4,:)-df(2,:).*df(3,:))./(df(4,:)-derprim.*df(3,:)));
	indf=find(~(dCs<0&dCs>df(1,:)&isfinite(dCs)));
	dCs(indf)=df(1,indf)/2;
	isoc.Cs(ind+1)=isoc.Cs(ind)+dCs;
	isoc.aws(ind+1)=isoc.aws(ind)+dCs.*derprim;
	%i;// a gauche des points de référence :
	ind=cumsum(nbs);
	derprim=1./betiso.derconc(isoc.aws(ind),betiso);
	df=[diff([isoc.Cs(ind+1);isoc.Cs(ind);isoc.Cs(ind-1)]);diff([isoc.aws(ind+1);isoc.aws(ind);isoc.aws(ind-1)])];
	dCs=real((df(1,:).*df(4,:)-df(2,:).*df(3,:))./(df(4,:)-derprim.*df(3,:)));
	indf=find(~(dCs<0&dCs>df(1,:)&isfinite(dCs)));
	dCs(indf)=df(1,indf)/2;
	isoc.Cs(ind)=isoc.Cs(ind+1)+dCs;
	isoc.aws(ind)=isoc.aws(ind+1)+dCs.*derprim;
	%i;// on coupe ce qui n'est pas utile.
	isoc.K=[diff(isoc.Cs)./(eps^5+diff(isoc.aws))];
	isoc.Cs=isoc.Cs(1:end-1);
	if abs(betiso.aw0flu-isoc.aws(1))>eps^.2;
	  error('isoc.aw0~=0 is not yet supported by sensc, bug in isothermeacoin or invisothermeacoin.')
	  %i;// useless workaround:
	  isoc.aws=isoc.aws-betiso.aw0flu;
	end
	isoc.aw0=0;
	return
      else
	error('a programmer bug64983.')
      end
    end
    derseconde=(betiso.conc(aw+eps.^.3,betiso)+betiso.conc(aw-eps.^.3,betiso)-2*betiso.conc(aw,betiso))/eps.^.6;
    %i;//	    betiso.awminawmax(1:2)=betiso.awminawmax(1:2)-eps.^.3;
    %i;//	    approxiso1=feval(param.iso.helperfunctions.approximate_isotherm,betiso);
    %i;//	    betiso.awminawmax(1:2)=betiso.awminawmax(1:2)+2*eps.^.3;
    %i;//	    approxiso2=feval(param.iso.helperfunctions.approximate_isotherm,betiso);
    %i;//	    betiso.awminawmax(1:2)=betiso.awminawmax(1:2)-eps.^.3;
    %i;// derseconde=(approxiso1.Cs+approxiso2.Cs-2*isoc.Cs)*eps^(-.6);
    aw=isoc.invconc(isoc.Cs,isoc);
    %i;// C0flu est tres mal gere, il faut mieux retoucher l'isotherme ...
    imin=max(find(aw<=betiso.aw0flu));
    if eps^.8+aw(imin)<betiso.aw0flu & aw(imin+1)>eps^.8+betiso.aw0flu
      isoc.Cs=[isoc.Cs(imin)+isoc.K(imin)*(betiso.aw0flu-aw(imin)),isoc.Cs(imin+1:end)];
      aw=[0,aw(imin+1:end)-betiso.aw0flu];
      derseconde=derseconde(imin:end);
    else
      isoc.Cs=[isoc.Cs(imin),isoc.Cs(imin+1:end)];
      aw=[0,aw(imin+1:end)-betiso.aw0flu];
      derseconde=derseconde(imin:end);
    end
    if ploton
      figure(82);
      subplot(222);
      hold on;
      plot(aw,isoc.Cs);
    end
    isoc.Cs=min([isoc.Cs(2:end),eps^(-7)],max([-eps^(-7),isoc.Cs(1:end-1)],isoc.Cs-(abs(nbinter(2))-1)*derseconde.*([0,diff(aw)]+[diff(aw),0]).^2));
    if ploton
      figure(82);
      subplot(223);
      plot(aw,derseconde);
      subplot(222);
      hold on;
      plot(aw,isoc.Cs,'x');
    end
    %i;// cummax.
    for i=2:length(isoc.Cs)
      isoc.Cs(i)=max(isoc.Cs(i),isoc.Cs(i-1)+eps^.9*(eps^5+abs(isoc.Cs(i-1))));
    end
    if ploton
      figure(82);
      subplot(222);
      hold on;
      plot(aw,isoc.Cs,'+');
    end
    isoc.K=[diff(isoc.Cs)./(eps^5+diff(aw)),isoc.K(end)];
    return
  end
  %i;// test of functions
  for aw=iso.awminawmax(1):diff(iso.awminawmax(1:2))/10:iso.awminawmax(2);
    if abs(iso.invconc(iso.conc(aw,iso),iso)-aw)>eps^.5
      disp([iso.conc(aw,iso),iso.invconc(iso.conc(aw,iso),iso),aw,eps^.5])
      error(['disagreement conv et invconc dans l''isotherme',iso.name])
    end
    if isfield(iso,'paramconc')
      if abs(iso.paramconc(iso.conc(aw,iso),aw,iso))>eps^.5
	error(['disagreement conc et paramconc dans l''isotherme',iso.name])
      end
    end
  end
  if iso.awminawmax(end)<-2;
    %i;// iso.awminawmax=[awmin,awmax,-2-1/numberofawintervals];
    iso.awminawmax=[iso.awminawmax(1:2),iso.awminawmax(1):-diff(iso.awminawmax(1:2))*(iso.awminawmax(end)+2):iso.awminawmax(2)];
  end
  isoc.Cs=iso.conc(iso.awminawmax,iso);
  if ploton
    disp(isoc.Cs)
  end
  %i;// quand on a défini diff par erreur, on peut debeuguer avec which diff, whos diff opuis finir avec clear diff.
  if iso.awminawmax(end)<0;
    %i;// iso.awminawmax=[awmin,awmax,-1/numberofCsintervals];
    isoc.Cs=isoc.Cs(1):-diff(isoc.Cs(1:2))*iso.awminawmax(end):isoc.Cs(2);
  end
  isoc.Cs=unique(isoc.Cs);
  isoc.aws=iso.invconc(isoc.Cs,iso);
  isoc.K=diff(isoc.Cs)./diff(isoc.aws);
  isoc.Cs=isoc.Cs(1:end-1);
  isoc.aw0=isoc.aws(1);
  if isoc.aw0>0
    %i;// patch ci-dessous, sinon ca ne marche pas, aw0 non nul est tres mal gere.
    isoc.aws=[0,isoc.aws];
    %i;// isoc.K=[isoc.Cs(1)/(min(isoc.Cs(1),isoc.aw0)),isoc.K];
    isoc.K=[.1,isoc.K];
    isoc.Cs=[isoc.Cs(1)-isoc.aw0*isoc.K(1),isoc.Cs];
    isoc.aw0=0;
  end
  %i;// fin du patch.
  isoc.approximates=iso;
  isoc.conc=isothermeacoin;
  isoc.invconc=invisothermeacoin;
%i;endfunction
function iso=baldev02STARCHandPE(name,partdestarchengparg,varargin)
  first=1;%i;first=[]
  if length(varargin)>0
    iso=varargin(1);
    iso=[iso{first}];
  end
  switch lower(name);%i;select lower(name);// compatibilité scilab matlab.
    case 'bet',
      iso.M_m=[4.54 3.35 2.52 1.70 0.62];
      iso.C=[12.51 16.87 16.23 18.30 28.85];
      iso.R2=[0.96 0.96 0.99 0.95 0.97];
      iso.awminawmax=[.1 .4];
    case 'smith',
      iso.M_b=[5.41 3.96 3.00 1.95 0.95];
      iso.M_a=[1.96 1.84 1.37 0.97 0.76];
      iso.R2=[0.96 0.98 0.97 0.99 0.99];
      iso.awminawmax=[0.3 0.9];
    case 'halsey',
      iso.a=[1.87 1.60 1.32 0.90 0.29];
      iso.b=[-0.19 -0.23 -0.22 -0.24 -0.31];
      iso.R2=[0.98 0.99 1.00 0.99 0.99];
      iso.awminawmax=[0.4 0.9];
    case 'caurie',
      iso.A=[4.47 3.25 2.48 1.64 0.81];
      iso.r=[-0.85 -0.98 -0.95 -0.97 -1.26];
      iso.R2=[0.99 0.99 0.97 0.96 0.9];
      iso.awminawmax=[0.3 0.9];
    case 'bradley',
      error('constants in the table II are inexact, please measure yourself the constants in figure 6 of this paper');
      %i;// NB je crois qu'il a oublié un log, et un changement de signe, d'après la figure 6.
      iso.K_1=exp([1.97 2.01 2.60 3.64 5.03]);
      iso.K_2=exp(100*[0.01 0.04 0.03 0.05 0.13]);
      iso.R2=[0.99 0.99 0.99 0.98 0.98];
      iso.awminawmax=[0.4 0.9];
    case 'oswin',
      iso.a=[7.04 5.39 4.06 2.65 1.46];
      iso.n=[0.15 0.18 0.18 0.20 0.28];
      iso.R2=[1.00 0.99 0.99 0.98 0.99];
      iso.awminawmax=[0.5 0.9];
  end
  f=fieldnames(iso);
  for param=f(1:3)';
    param=[param{first}];
    evalc(['iso.',param,'=iso.',param,'(round(6-partdestarchengparg*10));']);
  end
  iso.name=name;
%i;endfunction
function listest=allongelistest(listest,listechamp,resanalytic)
  first=1;%i;first=[]
  for st=listechamp
    for suffix={'','_energy'}
      evalc(strcat(['lengt=length(resanalytic.parametres',suffix{first},')']));
      for n=1:lengt
	listest(end+1)={strcat(['resanalytic.parametres',suffix{first},'(',num2str(n),').',st{first}])};
      end
    end
  end
%i;endfunction
function res=isothermeacoin
%i;//   C=iso.Cs(1)+sum(min([diff(iso.Cs),eps^(-2)]'*(1+0*aw),iso.K.*max(0,(1+0*iso.K)'*aw-iso.aw0-[0,cumsum(diff(iso.Cs)./iso.K(1:end-1))]'*(1+0*aw))));
res         = inline('iso.Cs(1)+sum(min([diff(iso.Cs),eps^(-2)]''*(1+0*aw),(iso.K''*(1+0*aw)).*max(0,(1+0*iso.K)''*aw-iso.aw0-[0,cumsum(diff(iso.Cs)./iso.K(1:end-1))]''*(1+0*aw))))','aw','iso');
%i;endfunction
function res=invisothermeacoin
res         = inline('iso.aw0+sum(min([diff(iso.Cs),eps^(-2)]''*(1+0*C),max(0,(1+0*iso.K)''*C-iso.Cs''*(1+0*C)))./(iso.K''*(1+0*C)))','C','iso'); % s = scale factor
  %i;// aw=iso.aw0+sum(min([diff(iso.Cs),eps^(-2)]'*(1+0*C),max(0,(1+0*iso.K)'*C-iso.Cs'*(1+0*C)))./(iso.K'*(1+0*C)))
%i;endfunction


%i;// TODO biblio a faire : 

%  google site:cat.inist.fr controlabilité estimation paramètres me donne :
%i;//   Modelling and control of an industrial PVC suspension polymerization reactor\n\nAuteur(s) / Author(s)\n\n   LEWIN D. R.^ (1) ;\n\nAffiliation(s) du ou des auteurs / Author(s) Affiliation(s)\n\n   ^(1) Department of Chemical Engineering, Technion, Haifa 32000, ISRAEL\n\nRésumé / Abstract\n\n   The modelling, parameter estimation and controllability analysis of an industrial batch PVC\n   suspension polymerization reactor is described. A mathematical model is developed in which the\n   parameters of a simplified reactor kinetics expression are fitted using plant data. The\n   resulting model can accurately predict the impact of initiator loading and coolant temperature\n   on the operability of the reactor.\n\nRevue / Journal Title\n\n   Computers & chemical engineering  (Comput. chem. eng.)  ISSN 0098-1354   CODEN CCENDW  Computers\n   and chemical engineering\n\nSource / Source\n\n   Congrès\n   European Symposium on Computer Aided Process Engineering_6. Part B\n   European Symposium on Computer Aided Process Engineering _ ESCAPE N^o6, Rhodes , GRECE\n   (26/05/1996)\n   Symposium of the Working Party on Computer Aided Process Engineering (CAPE) N^o29, Rhodes ,\n   GRECE (26/05/1996)\n   Event of the European Federation of Chemical Engineering (EFCE) N^o541, Rhodes , GRECE\n   (26/05/1996)\n   1996, vol. 20 (5 ref.), pp. S865-S870\n  ne fait pas de controle en ligne, seulement specification de parametres constants lors de l'expérience.
%i;// Remarks on uncertainty assessment and management in modeling and computation Auteur(s) / Author(s) BANKS H. T.^ (1) ; Affiliation(s) du ou des auteurs / Author(s) Affiliation(s) ^(1) Center for Research in Scientific Computation, Box 8205, North Carolina State University, Raleigh, NC 27695-8205, ETATS-UNIS Résumé / Abstract We discuss questions related to uncertainty in scientific computations for mathematical models.  A computationally tractable probabilistic framework to treat uncertainty in the estimation of parameters or inverse problems is given. The theory is illustrated by a simple computational example for the estimation of constant parameters in differential equations by treating the parameters as random variables.  Revue / Journal Title Mathematical and computer modelling  (Math. comput. model.)  ISSN 0895-7177 Source / Source Congrès Computation and Control VI Conference on Computation and Control N^o6, Bozeman, Montana , ETATS-UNIS (04/08/1998) 2001, vol. 33, n^o 1-3 (13 ref.), pp. 39-47
%i;// Nonlinearity, scale, and sensitivity for parameter estimation problems Auteur(s) / Author(s) GRIMSTAD Alv-Arne^ (1) ; MANNSETH Trond^ (2) ; Affiliation(s) du ou des auteurs / Author(s) Affiliation(s) ^(1) Department of Physics, University of Bergen, Allégt. 55, 5007 Bergen, NORVEGE ^(2) RF-Rogaland Research, Thormøhlensgt. 55, 5008 Bergen, NORVEGE Résumé / Abstract We study model nonlinearity and sensitivity for parameter estimation problems. Previous work has revealed a correlation between high nonlinearity, low sensitivity, and short scale for an ODE model. In this paper, we investigate whether this correlation is valid for a larger class of model functions. We set forth a proposition that says, in essence, that the correlation holds when the forward model output is an integral. Solutions to ODEs and PDEs can be viewed as integral models. If the proposition is true, it may explain an apparent conflict between earlier works on uncertainty analysis for parameter estimation problems. It could also have impact on the choice of solution algorithm for such problems. The validity of the proposition is assessed by studying nonlinearity and sensitivity for fairly general nonlinear models, and corresponding integrals of these models. Theory is developed and its predictions are tested through numerical examples. The focus is on general effects for the model classes (integrated/nonintegrated), i.e., effects that do not depend on any specific choice of model function within each class.  Revue / Journal Title SIAM journal on scientific computing  (SIAM j. sci. comput.)  ISSN 1064-8275   CODEN SJOCE3 Source / Source 2000, vol. 21, n^o6, pp. 2096-2113 (13 ref.) 
%i;// The problem of pole-zero cancellation in transfer function identification and application to adaptive stabilization Auteur(s) / Author(s) CAMPI M. C.^ (1) ; Affiliation(s) du ou des auteurs / Author(s) Affiliation(s) ^(1) Dipartimento di Elettronica per l'Automazione, Universita' di Brescia, via Branze, 38, 25123 Brescia, ITALIE Résumé / Abstract The asymptotic controllability of the identified system is a central problem in adaptive control. If controllability is ascertained, the analysis of even complex adaptive controllers based on multistep performance indices is drastically simplified. In this paper, we study the controllability issue in connection with the recursive least-squares (RLS) algorithm. We show that standard RLS does not generally provide models that are controllable. However, a variant of this method that preserves all the basic properties of the standard RLS and also guarantees asymptotic controllability is introduced. The algorithm can be safely used in any adaptive control system, provided that the control law is able to stabilize known invariant plants.  Revue / Journal Title Automatica  (Automatica)  ISSN 0005-1098   CODEN ATCAA9
%i;// Modelling extrusion cooking Auteur(s) / Author(s) LI Chin-Hsien^ (1) ; Affiliation(s) du ou des auteurs / Author(s) Affiliation(s) ^(1) CSIRO Mathematical and Information Sciences, Locked Bag 17, North Ryde, NSW 2113, AUSTRALIE Résumé / Abstract A one-dimensional computer model of extrusion cooking has been developed. This model can simulate and predict extruder behaviour (such as pressure, temperature, fill factor, residence time distribution, shaft power, degree of cook) under various operating conditions (such as feed rate, screw speed, feed temperature/moisture, barrel temperature). With a very fast and efficient solution algorithm, this model runs fast on a PC, taking only a fraction of a second to a few seconds. In addition, with a graphic user interface, this model is user friendly, and can be used as a tool for designing, optimizing, and controlling extruder operating conditions.  Some results of numerical simulation will be presented.  Revue / Journal Title Mathematical and computer modelling  (Math. comput. model.)  ISSN 0895-7177 Source / Source Integrated modelling and estimation 2001, vol. 33, n^o 6-7 (8 ref.), pp. 553-563 
%i;//  Direct and inverse inequalities for the isotropic Lamé system with variable coefficients Auteur(s) / Author(s) GRASSELLI M.^ (1) ; IKEHATA M.^ (2) ; C M. Y.^ (3) ; Affiliation(s) du ou des auteurs / Author(s) Affiliation(s) ^(1) Dipartimento di Matematica F. Brioschi, Politecnico di Milano, Via E. Bonardi 9, 20133 Milano, ITALIE ^(2) Department of Mathematics, Faculty of Engineering, Gunma University, Kiryu 376-8515, JAPON ^(3) Department of Mathematical Sciences, The University of Tokyo, 3-8-1 Komaba Meguro, Tokyo 153, JAPON Résumé / Abstract Considérons un problème de Cauchy-Dirichlet pour le système de Lamé aux coefficients variables dans le cas isotrope. Nous déterminons une estimation pour la norme L[2] de la traction surfacique par rapport à les données initiales et à la force de volume. Nous demontrons par la suite que l'energie élastique, en absence des forces de volume, peut être contrôlée par la norme L[2] de la traction surfacique exercée sur une partie convenable du bord, à condition que l'intervalle de temps soit suffisament grand. Ces inégalités sont fondamentales pour l'application de la méthode nommée HUM (Hilbert Uniqueness Method). Elles peuvent aussi être utilisées pour résoudre un problème inverse pour le système de Lamé.  Revue / Journal Title Comptes rendus de l'Académie des sciences. Série IIb, Mécanique  (C. r. Acad. sci., Sér. IIb Méc.)  ISSN 1620-7742 
%i;// Adaptive pole positioning in MIMO linear systems by periodic multirate-input controllers Auteur(s) / Author(s) ARVANITIS K. G.^ (1) ; KALOGEROPOULOS G.^ (2) ; SANTAS E. A.^ (3) ; Affiliation(s) du ou des auteurs / Author(s) Affiliation(s) ^(1) Department of Electrical and Computer Engineering, National Technical University of Athens, Division of Computer Science, Zographou 15773, Athens, GRECE ^(2) Department of Mathematics, Section of Mathematical Analysis, University of Athens, Panepistimiopolis 15784, Athens, GRECE ^(3) Ministry of Interior, Administration and Decentralization, Region of Attica, 60 Theras Str., 11262 Athens, GRECE Résumé / Abstract In this paper, the certainty equivalence principle is used to combine the identification method with a control structure derived from the pole placement problem, which rely on periodic multirate-input controllers. The proposed adaptive pole placers, contain a sampling mechanism with different sampling period to each system input and rely on a periodically varying controller which suitably modulates the sampled outputs and reference signals of the plant under control. Such a control strategy allows us to arbitrarily assign the poles of the sampled closed-loop system in desired locations and does not make assumptions on the plant other than controllability and observability of the continuous and the sampled system, and the knowledge of a set of structural indices, namely the locally minimum controllability indices of the continuous-time plant. An indirect adaptive control scheme is derived, which estimates the unknown plant parameters (and consequently the controller parameters) on-line, from sequential data of the input and outputs of the plant, which are recursively updated within the time limit imposed by a fundamental sampling period T[0]. Using the proposed algorithm, the controller determination is based on the transformation of the discrete analogous of the.  Revue / Journal Title Journal of mathematical analysis and applications  (J. math. anal. appl.)  ISSN 0022-247X CODEN JMANAK Source / Source 1999, vol. 237, n^o2, pp. 464-504 (40 ref.) 
%i;// TODO apprendre a appeler matlab via wine avec /C/MATLAB6p5/extern/examples/refbook/ ...
%i;//  vim:softtabstop=2:shiftwidth=2:guioptions+=a
%i;//  :noremap à :execute 'w'<bar>:w! Y:\common\senscrep20080616\<c-r>%
