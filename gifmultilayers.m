function out = gifmultilayers(varargin)
%MULTLAYERS generate GIF movies for multilayer structures (default = ABC + Food)
% SYNTAX: gifmultilayers('property1',value1,'property2',value2,...,'keyword1','keyword2',...)
%         out = gifmultilayers(...)
%
% LIST OF PROPERTY/VALUE: see code and default object (not fully documented at this stage)
%       default properties can be retrieved with: S = gifmultilayers('initialize')
%       property inheritance from Sold can be obtained by placing Sold at the last argument to be overdefined as:
%       S = gifmultilayers('initialize','propertytobeoverdefined',value....,Sold)
%
% LIST OF KEYWORDS
%           'initialize'    returns inputs
%           'restart'       overwrites any prefetch file
%           'ploton'        plots frames
%           'makegif'       generate a GIF file from plotted frames
%           
% An early version of this code can be used in: 
%   "C:\Data\Olivier\INRA\Etudiants & visiteurs\Mai Nguyen\soutenance\base\desorption\desorption.m"
%
% OUTPUTS (not fully documented at this stage)
%
%
% ========================================================================================================
% EXAMPLES to illustrate the principles of transfer between multilayers (generated for Saint-Gobain USA)
% ========================================================================================================
%{
% ======= Reference ======= 
gifmultilayers('makegif')
% ==== functional barrier =====
gifmultilayers('D',[1 1 .1 1e4],'filename','funcbarrier','makegif')
% ==== functional barrier =====
gifmultilayers('D',[1 1 .1 1e4],'k',[1 1 50 1],'Fo',logspace(-2,log10(200),1200),'filename','funcbarrier2','makegif')
% ==== functional barrier =====
gifmultilayers('D',[1 1 .1 1e4],'k',[1 1 50 20],'Fo',logspace(-2,log10(200),1200),'filename','funcbarrier3','makegif')
% ==== many layers =====
manymodel = struct(...
    'n', [1.4  0.6 1  .5  1/20]*20,...
    'l', [1    1   1   1  100],...  
    'C0',[.2   0   1   0  0],...
    'D', [2    1   1  .2  1e4 ],...
    'k', [3   1   1/3  1   1/3] );
position = [200   150   1280   800];
gifmultilayers(manymodel,'Fo',logspace(-2,log10(200),1200),'filename','manylayers','makegif','position',position,'markersize',4,'title','\bfQUADRIlayer ABCD\rm')
%}
%
% ====================================================================================================
% EXAMPLES to illustrate superposition effects (generated for Sartorius-Stedim, France, Germany)
% ====================================================================================================
%{
    % reference [F1 P F2], note the asymmetry
    local = pwd; % set your path here
    common = struct('nmesh',400,'Fo',linspace(0,sqrt(25),1e3)'.^2,'local',local); % Fo based on the smallest element (not layer), Fo for layer 2=25/20=1.25
    Sref = gifmultilayers('initialize','filename','superposition_reference',...
                          'n',[1 20 1],'l',[20 1 30],'C0',[0 1000 0],'D',[1e4 1 1e4],'k',[3 1 2],common);
    Rref = gifmultilayers(Sref);
    figure, hp = plot(Rref.Fo_input,Rref.Cmean); xlabel('Fo'), ylabel('Conc.'); legend(hp,{'F_1' 'P' 'F_2'})
    hold on, plot(xlim,Rref.Ceq([1 1],:),':','linewidth',2)
    % left F (F1) only
    Sleft = gifmultilayers('initialize','restriction',1:2,'filename','superposition_left',Sref); % matches [F1 P]
    Sleft.C0(2) = Sleft.C0(2) * Rref.Ceq(1)/Sleft.Ceq(1); % update the source (position 1)
    Rleft =  gifmultilayers(Sleft);
    % right F (F2) only
    Sright = gifmultilayers('initialize','restriction',2:3,'filename','superposition_right',Sref); % matches [P F2]
    Sright.C0(1) = Sright.C0(1) * Rref.Ceq(3)/Sright.Ceq(2); % update the source (position 2)
    Rright =  gifmultilayers(Sright);
    % plot comparison
    figure, hold on, hp = [
            plot(Rref.Fo_input,Rref.Cmean(:,1),'-','color',rgb('Teal'),'linewidth',2)
            plot(Rref.Fo_input,Rref.Cmean(:,3),'-','color',rgb('Crimson'),'linewidth',2)
            plot(Rleft.Fo_input,Rleft.Cmean(:,1),'-','color','b','linewidth',1)
            plot(Rright.Fo_input,Rright.Cmean(:,2),'-','color','r','linewidth',1)
            plot(Rref.Fo_input,Rref.Cmean(:,2),'-','color','k','linewidth',2)
            plot(Rref. Fo_input, (interp1(Rleft.Fo_input,Rleft.Cmean(:,2),Rref.Fo_input) + ...
                                 interp1(Rright.Fo_input,Rright.Cmean(:,1),Rref.Fo_input))/2,...
              '--','color','k','linewidth',1)
                ];
    xlabel('Fo'), ylabel('Conc.'); legend(hp,{'exact F_1' 'exact F_2' 'F_1 alone' 'F_2 alone' 'exact P' 'average of 2 approximations'})
    % Concentration profile for a particular Fo
    Fotarget =2; it = nearestpoint(Fotarget,Rref.Fo_input);
    figure 
    subplot(131), hold on, Rref.plotbnd(),Rref.plotpart(it), Rref.plot(it,it), axis tight, title('\bfexact solution')
    subplot(132), hold on, Rleft.plotbnd(),Rleft.plotpart(it), Rleft.plot(it,it), axis tight, title('left approximation')
    subplot(133), hold on, Rright.plotbnd(),Rright.plotpart(it), Rright.plot(it,it), axis tight, title('right approximation')
%}
%
% As above but generate GIF images (please run the example above before launching this code)
%{
    gifmultilayers('xticklabel',{'0' '1' '2' '3'},'title','exact solution',Sref,'makegif')
    gifmultilayers('xticklabel',{'0' '1' '2'},'title','left approximation',Sleft,'makegif')
    gifmultilayers('xticklabel',{'1' '2' '3'},'title','right approximation',Sright,'makegif')
%}
%
% ======================================================================
% alternative scenario based on mass balance and thickness correction
% Note: that additional sources and reflections are required to make the
%       approximation exact (see applications of Green's transform)
% =======================================================================
%{
    meq = Rref.Ceq .* Sref.l;
    ratio = meq([1 3])./sum(meq([1 3])); % main assumption: each side of P contributes proportionally to the final eq value
    S2left = gifmultilayers('initialize','restriction',1:2,'filename','superposition2_left',Sref); % matches [F1 P]
    S2left.l(2) = S2left.l(2) *ratio(1); % update the source thickness (position 1)
    R2left =  gifmultilayers(S2left);
    S2right = gifmultilayers('initialize','restriction',2:3,'filename','superposition2_right',Sref); % matches [P F2]
    S2right.l(1) = S2right.l(1) *ratio(2); % update the source thickness (position 1)
    R2right =  gifmultilayers(S2right);
    figure, hold on, hp = [
        plot(Rref.Fo_input,Rref.Cmean(:,1),'-','color',rgb('Teal'),'linewidth',2)
        plot(Rref.Fo_input,Rref.Cmean(:,3),'-','color',rgb('Crimson'),'linewidth',2)
        plot(Rleft.Fo_input,R2left.Cmean(:,1),'-','color','b','linewidth',1)
        plot(Rright.Fo_input,R2right.Cmean(:,2),'-','color','r','linewidth',1)
        plot(Rref.Fo_input,Rref.Cmean(:,2),'-','color','k','linewidth',2)
        plot(Rref. Fo_input, ratio(1)*interp1(R2left.Fo_input,R2left.Cmean(:,2),Rref.Fo_input) + ...
                              ratio(2)*interp1(R2right.Fo_input,R2right.Cmean(:,1),Rref.Fo_input),...
        '--','color','k','linewidth',1)
        ]; title('2nd strategy: mass balance + thickness correction'), xlabel('Fo'), ylabel('Conc.'); 
    legend(hp,{'exact F_1' 'exact F_2' 'F_1 alone' 'F_2 alone' 'exact P' 'average of 2 approximations'})
    figure, Fotarget =2; it = nearestpoint(Fotarget,Rref.Fo_input); % Concentration profile for a particular Fo
    subplot(131), hold on, Rref.plotbnd(),Rref.plotpart(it), Rref.plot(it,it), axis tight, title('\bfexact solution')
    subplot(132), hold on, R2left.plotbnd(),R2left.plotpart(it), R2left.plot(it,it), axis tight, title('left approximation')
    subplot(133), hold on, R2right.plotbnd(),R2right.plotpart(it), R2right.plot(it,it), axis tight, title('right approximation')
%}
% As above but generate GIF images (please run the example above before launching this code)
%{
    gifmultilayers('xticklabel',{'0' '1' sprintf('%0.3g',1+ratio(1)) },'title','left approximation (2)',S2left,'makegif')
    gifmultilayers('xticklabel',{sprintf('%0.3g',1+ratio(1)) '2' '3'},'title','right approximation (2)',S2right,'makegif')
%}

% Migration 2.1 - 17/12/2015 - INRA\Olivier Vitrac - rev. 21/04/2017

% Revision history
% 17/12/2015 release candidate as a function with examples, based on the script used by Mai and a former script used for EU Commission in 2012
% 23/12/2015 consolidation and output, documentation is still succinct but sufficient for advanced users
% 24/12/2015 function renamed as gifmultilayers (to be integrated within MIGRATION toolbox)
% 24/12/2015 example for Sartorius completed, add restriction, plots with arbitrary concentration units
% 25/12/2015 second approximation of superpoposition (preferred): mass balance and thickness correction
% 26/02/2017 returns npart() for 3D
% 21/04/2017 add nscale and nscaleleg to remove the effect of n (legacy maintained)

% default
ndefault = 20;
nFodefault = 400;
Fomin = 1e-2; Fomax = 200;
default = struct(...
    'local',pwd,...
    'filename','multilayer',...
    'N0',100,...
    'C0',[0 1 0 0],... ABC + Food
    'l',[1 1 1 100],... ABC + Food
    'D',[1 1 1 1e4],... ABC + Food
    'k',[1 1 1 1],... ABC + Food
    'n',[ndefault ndefault ndefault 1],... ABC + Food
    'nscale',1,... 1 is for legacy
    'nscaleleg',2,...
    'Fo',logspace(log10(Fomin),log10(Fomax),nFodefault),...
    'ind',[],...
    'colors', rgb({'Salmon' 'FireBrick'}),...
    'colormap',@jet,...
    'colormapbnd',@summer,...
    'delaytime',1/20,...
    'position',[222   261   800   600],...
    'title','\bfTRIlayer ABC\rm',...
    'linewidth',3,...
    'markersize',2,...
    'fontsize',14,...
    'xticklabel','',...
    'nmesh',1200,...
    'restriction',[] ...
    );
keywords = {'makegif' 'restart' 'ploton' 'initialize'};

%% arg check
o = argcheck(varargin,default,keywords,'case');
prefetchfile = fullfile(o.local,sprintf('%s.mat',o.filename));
giffile = fullfile(o.local,sprintf('%s.gif',o.filename));
o.n = o.n(:)';
o.l = o.l(:)';
o.D = o.D(:)';
o.k = o.k(:)';
o.C0 = o.C0(:)';
if isempty(o.n), o.n=ndefault; end
if length(o.n)==1, o.n = ones(size(o.l))*o.n; end
nlayers = unique([length(o.n) length(o.l) length(o.D) length(o.k) length(o.C0)]);
if length(nlayers)>1, error('the size of n,l,D,k,C0 are inconsistent (min=%d max=%d)',min(nlayers),max(nlayers)), end
if ~isempty(o.restriction)
    if any(setdiff(o.restriction,1:nlayers))
        error('The property ''restriction'' is applied outside existing layers 1..%d', nlayers)
    end
    o.n = o.n(o.restriction);
    o.l = o.l(o.restriction);
    o.D = o.D(o.restriction);
    o.k = o.k(o.restriction);
    o.C0 = o.C0(o.restriction);
end
if o.makegif || ~nargout, o.ploton = true; end
if isempty(o.xticklabel)
    o.xticklabel = [arrayfun(@(k) sprintf('%d',k),0:length(o.k)-1,'UniformOutput',false) {'food'}];
end
if o.initialize, out=rmfield(o,keywords); out.restriction=[]; out.Ceq = eqConc(out); return, end

%% functions
markersize = o.markersize;
N0 = o.N0;
colbasic = o.colors;
colorbnd = o.colormapbnd;
hcolor   = @(h,col) arrayfun(@(i)set(h(i),'color',col(i,:),'linestyle','-','linewidth',1.2),1:length(h));
npart    = @(C,dx) round(N0*C*dx);    % number of particles
x        = @(s) cumsum([0 s(2,:)],2); % positions
plotbnd  = @(s,Cmax)  hcolor(line(repmat(x(s),2,1),[zeros(1,size(s,2)+1);Cmax*ones(1,size(s,2)+1)]),colorbnd(size(s,2)+1));
plotpart = @(C,x,it,colidx,Cmax) arrayfun(@(C,x,dx,icol) ... Cobject positions it colorindices
                plot(x+rand(npart(C,dx),1)*dx,Cmax*rand(npart(C,dx),1),'o','color',colbasic(icol,:),'markerfacecolor',colbasic(icol,:),'markersize',markersize),...
                C.C(it,:)/Cmax,x(1:end-1),diff(x),colidx);
%% SIMULATION
if ~exist(prefetchfile,'file') || o.restart
    s = [
        expandmat(o.C0,o.n) %ones(1,n)   ones(1,n)    zeros(1,n) 0
        expandmat(o.l,o.n) %ones(1,n)    ones(1,n)    ones(1,n)  100
        expandmat(o.D,o.n) %ones(1,n) 10*ones(1,n)    ones(1,n)  1e4
        expandmat(o.k,o.n) %ones(1,n)    ones(1,n)    ones(1,n)  1
        ];
    n = max(o.n);
    Fo = o.Fo*n.^o.nscale;
    dispf('run simulation, defined as:'),disp(s)
    C = roe_patankar(s,Fo,0,'','struct',1,'nmesh',o.nmesh);
    save(prefetchfile,'o','s','Fo','C','n')
    dispf('generate the following prefetch file:')
    fileinfo(prefetchfile)
else
    dispf('use the following prefetch file:')
    fileinfo(prefetchfile)
    otmp = o;
    load(prefetchfile)
    o = argcheck(otmp,o);
end

%% plot profiles and animate
nFo = length(Fo);
if isempty(o.ind), itlist = 1:nFo; else, itlist = o.ind((o.ind>=1) & (o.ind<=nFo)); itlist = round(itlist(:)'); end
colfull = o.colormap(512); colfull = colfull(round(interp1(linspace(0,sqrt(Fo(end)-Fo(1)),512).^2,1:512,Fo-Fo(1),'pchip')),:);
col = colfull(itlist,:); j = 1;
if o.ploton
    Cmax = max(o.C0);
    formatfig(figure,'figname','profiles','paperposition',[4.3270    9.2298   12.3301   11.2178],'color','w','position',o.position)
    drawnow
    for it = itlist
        clf, drawnow, hold on
        plotbnd(s,Cmax), drawnow
        plotpart(C,x(s),it,(s(1,:)>0)+1,Cmax), drawnow
        plot(C.x,C.Cx(it,:),'-','linewidth',o.linewidth,'color',col(j,:)); drawnow
        ylabel('concentration profile','fontsize',o.fontsize) % (C_F^{t\rightarrow\infty}/C_P^{t=0}\approx1/5)
        xlabel('dimensionless position x/l_0','fontsize',o.fontsize)
        formatax(gca,'fontsize',13,'ylim',[0 Cmax],...
            'xlim',[0 sum(s(2,:))],...
            'xtick',[0 cumsum(o.n.*o.l)],...
            'xticklabel',o.xticklabel)
        title(sprintf('%s   Fo = D\\cdott/l_0^2 = \\bf%s\\rm \\fontsize{10}(-)',o.title, formatsci(Fo(it)/n^o.nscaleleg)),'fontsize',o.fontsize)
        drawnow
        % make gif
        if o.makegif
            frame = getframe(1); im = frame2im(frame); [A,map] = rgb2ind(im,256);
            if it == itlist(1)
                imwrite(A,map,giffile,'gif','LoopCount',Inf,'DelayTime',o.delaytime);
            else
                imwrite(A,map,giffile,'gif','WriteMode','append','DelayTime',o.delaytime);
            end
        end
    end
end

%% output
if nargout
    out = C;
    out.prefetchfile = prefetchfile;
    out.ndivisions = n;
    out.Fo_input = out.Fo/n;
    out.ind = arrayfun(@(start,stop) start:stop,cumsum([1 o.n(1:end-1)]),cumsum(o.n),'UniformOutput',false);
    tmp = cellfun(@(j) mean(C.C(:,j),2),out.ind,'UniformOutput',false);
    out.Cmean = cat(2,tmp{:});
    out.Ceq   = eqConc(o);
    Cmax = max(o.C0);
    out.npart = @(N0,C,dx) round(N0*C*dx);
    out.plotbnd = @() plotbnd(s,Cmax);
    out.plotpart = @(it) plotpart(C,x(s),it,(s(1,:)>0)+1,Cmax);
    out.plot = @(it,icolor) plot(C.x,C.Cx(it,:),'-','linewidth',o.linewidth,'color',col(icolor,:));
end

end % end function

%% PRIVATE FUNCTION
function Ceq = eqConc(p)
    % equilibrum concentration for the last layer (as if it was F)
    l_n = p.l.*p.n;
    CFeq = sum( (l_n/l_n(end)) .* p.C0 ) ./ sum( (p.k(end)./p.k) .* (l_n/l_n(end)) );
    Ceq   = CFeq * p.k(end)./p.k;
end
