function [Cout,Fout,Cfullout,resout]=roe_patankar(stack,Fo,ploton,printname,varargin)
%ROE_PATANKAR simulates a generalized Roe test using the patankar method (see SENSPATANKAR for help)
%   syntax: [Cout,Fout,Cfullout,s]=roe_patankar([stack],[Fo],[ploton],[printname],['keyword','value'])
%   inputs:
%       stack     = stack pattern (1..4)xnlayer array
%                   with nlayer = number of layers
%                   row 1 (mandatory) = concentrations in each layer (usually a combination of 0 and 1)
%                   row 2 = dimensionless thickness of each layer (usually filled of 1)
%                   row 3 = dimensionless diffusion coefficient (usually filled of 1)
%                   row 4 = dimensionless activity/Henry coefficient (usually filled of 1)
%       Fo        = dimensionless Fourrier time (respectively to a layer of thickness 1 with a diffusion coefficient of 1
%                   (default values = 0:.1:10)
%       ploton    = plot flag (default value = 0)
%       printname = name of the figure to be printed as a pdf file. No printing if empty (default)
%       recognized keywords: "nmesh' (default=300), 'options' (default, see code), 'struct' (if true output as a structure)
%   advanced inputs: 'iref' (see senspatankar) and 'eval' (see applesim)
%   outputs:
%           Cout = nFoxnlayer numeric array, concentration in each layer (averaged over the thickness)
%           Fout = 1xnFo numeric array, Fourier times where solutions are available
%        Cfullout= structure with fields
%           Cfullout.x = nnodesx1 discretization nodes of the finite volume method
%           Cfullout.C = nFoxnnodes "continuous" concentrations profiles
%             res= ouput of senspatankar
%
% Example:
%   roe_patankar([],[],1,'roe_patankar_maxorder_5')
                                                                                                                 
% MOISAN TOOLBOX 1.0 - 12/03/09 - Guillaume Gillet - Olivier Vitrac - rev. 19/01/2019

% revision history
% 16/06/09 tunning
% 25/05/10 remove simulation from loop
% 03/06/10 add variable D and k between layers
% 06/09/10 add property nmesh, options
% 13/02/11 add 'struct' as plot
% 13/11/11 add resout
% 14/06/12 accept Fo values beyond 1000
% 25/07/13 return respatankar.F
% 05/12/17 add iref, eval
% 16/01/19 fix argcheck for options
% 19/01/19 add L, robusteps

% definitions
stack_default = [0 0 0 0 1 0 0 0 0 1 0 0 0 0];
Fo_default = [.01:.01:.05 .08 .1:.1:1 1.5 2 5 10];
ploton_default = false;
options_default = odeset('RelTol',1e-4,'AbsTol',1e-4,'Initialstep',1e-8,'Maxstep',.1,'Maxorder',2); %options of ode used in SENSPATANKAR
prop_default = struct('nmesh', 300,'struct',false,'iref',[],'eval',[],'L',1);
% options		= odeset('RelTol',1e-4,'AbsTol',1e-4,'Initialstep',1e-8,'Maxstep',.01,'Maxorder',5); % options of ode used in SENSPATANKAR
robusteps = 100*eps;

% arg check
%if nargchk(0,4,nargin), error('syntax: [Cout,Fout,Cfullout]=roe_patankar([stack],[Fo],[ploton],[printname])'), end %#ok<NCHK>
if nargin<1, stack = []; end
if nargin<2, Fo = []; end
if nargin<3, ploton = []; end
if nargin<4, printname = []; end
[options,remain] = argcheck(varargin,struct('options',options_default),'','nostructexpand','keep');
prop = argcheck(remain,prop_default);
prop.options = options.options;
if isempty(stack), stack = stack_default; end
if ndims(stack)==1, stack = stack(:)'; end
nlayer = size(stack,2);
if size(stack,1)<2, stack(end+1,:) = ones(1,nlayer); end
if size(stack,1)<3, stack(end+1,:) = ones(1,nlayer); end
if size(stack,1)<4, stack(end+1,:) = ones(1,nlayer); end
if isempty(Fo), Fo = Fo_default; end
if isempty(ploton), ploton = ploton_default; end

% From Guillaume (2009) // REMOVED by OV June 2010
% nstacknorm = sum(stack(2,:)); % total thickness
% C0 = NaN(1,nstacknorm);
% iC = 1;
% for nstack=1:nstacknorm % conversion in layers of thickness equal to 1
%     if nstack<=sum(stack(2,1:iC))
%         C0(nstack) = stack(1,iC);
%     else
%         iC = iC+1;
%         C0(nstack) = stack(1,iC);
%     end
% end
% // end REMOVED section

% Additional parameters ()


% problem formulation
F = struct('Bi'		 , 0,...	Biot [hm.L1/D]
           'k'		 , stack(4,:),...[0.5 3 2],...	ki, i=1 (layer in contact with the liquid)
           'D'       , stack(3,:),... diffusion coefficient
           'k0'      , 1,... 0 = liquid
           'l'       , stack(2,:),...[50 20 10 120]*1e-6,... m
		   'L'		 , prop.L,...	dilution factor (respectively to iref)
		   'C0'		 , stack(1,:),...	initial concentration in each layer
           'nmesh'   , prop.nmesh,...
		   'options' , prop.options);
if ~isempty(prop.iref), F.iref = prop.iref; end

% Simulation
nFo = length(Fo);
Fomax = max(Fo);
if Fomax<=0.01,  F.t = 0:.00001:Fomax;
elseif Fomax>1000, F.t = [0:.0001:.05 .1:.1:.9 1:1:9 10:10:100 linspace(200,Fomax,100)];
elseif Fomax>100, F.t = [0:.00001:.005 .01:.01:10 10.1:.1:10 10.5:.5:20 21:1:50 52:2:100 100:5:Fomax];
elseif Fomax>50, F.t = [0:.00001:.005 .01:.01:10 10.1:.1:10 10.5:.5:20 21:1:50 52:2:Fomax];
elseif Fomax>20, F.t = [0:.00001:.005 .01:.01:10 10.1:.1:10 10.5:.5:Fomax];
elseif Fomax>10, F.t = [0:.00001:.005 .01:.01:10 10.1:.1:Fomax];
else F.t = [0:.00001:.005 .01:.01:Fomax];
end
F.t = unique([Fo(:)' F.t]);
respatankar = senspatankar(F);
% time rescaling
respatankar.t = respatankar.t * respatankar.timebase;
% geometric interpretation
xlayerpos = [0 cumsum(F.l)];
xmesh = xlayerpos(end) * respatankar.x/respatankar.x(end); % length rescaling
indlayer = false(length(xmesh),nlayer);
left = -Inf;
for ilayer=1:nlayer
    indlayer(:,ilayer) = (xmesh>left) & (xmesh<=(xlayerpos(ilayer+1)+robusteps));
    left = max(xmesh(indlayer(:,ilayer)));
end
% averaging for each layer
res = repmat(struct('C',[],'Cx',[],'layer',[]),nFo,1);
for i=1:nFo
    iFo = (F.t==Fo(i));
    res(i).C = NaN(1,nlayer);
    res(i).Cx = respatankar.Cx(iFo,:);
    for ilayer=1:nlayer
        res(i).C(ilayer) = trapz(xmesh(indlayer(:,ilayer)),respatankar.Cx(iFo,indlayer(:,ilayer)))/F.l(ilayer);
    end
end

% figure (if asked)
if ploton
    figure('paperposition',[0.2 0.184 27.8 20.6],'paperorientation','landscape');
    hs = subplots(ones(1,ceil(sqrt(nFo))),ones(1,ceil(nFo/ceil(sqrt(nFo)))),.06,.08,'alive',1:nFo)';
    for i=1:nFo
        subplot(hs(i)), hold on
        bar(xlayerpos(1:end-1)+F.l/2,cat(1,res(i).C))
        plot(xmesh,res(i).Cx)
        set(gca,'ylim',[0 1.1*max([res(i).C])],...
                'xlim',xlayerpos([1 end]),...
                'xtick',.5:5:15.5,...
                'xticklabel',{'0' '5' '10' '15'})
        ht = title(sprintf('Fo = %g',Fo(i))); set(ht,'fontsize',11)
        xlabel('x','fontsize',10)
    end
end

% print (if asked)
if ~isempty(printname), print_pdf([],printname,cd), end

% nargout
C     = cat(1,res.C);
Cfull = cat(1,res.Cx);
if nargout, Cout = C; end
if nargout>1, Fout = Fo; end
if nargout>2
    Cfullout.x = xmesh;
    Cfullout.C = Cfull;
end
if prop.struct
    Cout = struct('Fo',Fo,'C',C,'x',xmesh,'Cx',Cfull,'F',respatankar.F);
end
if nargout>3, resout = respatankar; end

% eval if required
if ~isempty(prop.eval)
    if isa(prop.eval,'function_handle')
        Cout = prop.eval(Cout);
        return
    else
        error('eval must define an anonymous function')
    end
end