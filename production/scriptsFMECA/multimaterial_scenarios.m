%% SCRIPT TO ASSESS MODAL CRITICALITY of multimaterial scenarios

% INRA\Olivier Vitrac - 29/08/11    rev:
% Revision history :
% 30/08/11 "sim" multimaterial_scenarios_V3
% 31/08/11 add interpretation
% 31/08/11 "sim" multimaterial_scenarios_V4 (change in affinity of migrants)
% 31/08/11 "sim" multimaterial_scenarios_V5 (change in concentration of X, time of hotfilling step)
% 31/08/11 "sim" multimaterial_scenarios_V6
% 31/08/11 add 4x2 pareto charts

clear all

%% Definition
switch localname
    case 'WSLP-OLIVIER2'
        local = 'D:\Data\Olivier\INRA\Etudiants & visiteurs\Audrey Goujon\production\';
%         physicochemfolder = fullfile(find_path_toolbox('migration'),'database');
    case 'WS-MOL4'
        local = 'C:\Data\Olivier\Audrey_Goujon\Matlab_FMECA\production';
%         physicochemfolder = fullfile(find_path_toolbox('migration'),'database'); %can be setup to another location
    case 'mol10.agroparistech.fr'
        local = '/home/olivier/MaiNguyen/ANR_SFPD/WP2_FMECA/production';
%         physicochemfolder = fullfile(find_path_toolbox('migration'),'database');        
end

% Inputs (4 simulations)
nsim = 4; % number of simulations
scenariofolder = fullfile(local,'inputs');
fmecamainfile = 'multimaterial_scenarios_V6.ods'; % file containing sim table
fmecasheetname = arrayfun(@(isim) sprintf('sim%d',isim),1:nsim,'UniformOutput',false); % assuming sim1, sim2, sim3,...

% Paths for results and modal database (output)
% references simulations stored in (relatively to local):  'multilayer_scenarios\reference1'    'multilayer_scenarios\reference2'  ...
outputpathref = arrayfun(@(isim) fullfile('multilayer_scenarios',sprintf('reference%d',isim)),1:nsim,'UniformOutput',false);  % folder that contains all reference outputs
% with deletions (relatively to local): 'multilayer_scenarios\omission1'    'multilayer_scenarios\omission2 ....
outputpathdel = arrayfun(@(isim) fullfile('multilayer_scenarios',sprintf('omission%d',isim)),1:nsim,'UniformOutput',false); % folder that contains all omission outputs
figurefolder = fullfile(local,'figures'); if ~exist(figurefolder,'dir'), mkdir(figurefolder); end
fmecadbfile = 'fmecamodal.mat'; % base of conditions

% Solver parameters
nmesh = 200;   % refined 500
RelTol = 1e-4; % refined 1e-6
AbsTol = 1e-4; % refined 1e-6
options = odeset('RelTol',RelTol,'AbsTol',AbsTol,'Initialstep',1e-8,'Maxstep',1e3,'Maxorder',2);

%% simulations (restart not required for interpretation)
for isim=1:nsim
    % Reference simulations
    fmecadef = struct( ...
        'local',local,...
        'fmecasheetname',fmecasheetname{isim},...
        'fmecamainfile',fmecamainfile,... FMECAengine assumes always %local%/inputs
        'database',struct([]),... dummy database (to enable the use of a key such as Dpiringer('PP',182,25)')
        'inputpath','inputs',...
        'outputpath',outputpathref{isim}, ...
        'fmecadbfile',fmecadbfile,...
        'nmesh',nmesh,...
        'options',options ...
        );
    [fmecadbref,rawdata,dataref] = fmecaengine(fmecadef); % dataref (with inherited fields) is keep for modifications
    % Omission
    data = dataref;
    [data(1:length(data)).Foscale] = deal('[0 0]'); % add deletion "code" to all rows
    fmecadb=fmecaengine('local',local,...
        'fmecasheetname',fmecasheetname{isim},...
        'fmecamainfile',data,...
        'inputpath','omission.ods',...
        'outputpath',outputpathdel{isim}, ...
        'fmecadbfile',fmecadbfile, ...
        'nmesh',nmesh,...
        'options',options ...
        );

end

%% interpretation (based on CF)
figname = [regexprep(fmecamainfile,'\..*$','') '_reference']; % figure name
severity=@(CF,SML) 100*99./max(100*SML./CF-1,0); % Severity definition
mollist = {'X' 'Y'}; % list all molecules
nmol = length(mollist);
col = rgb({'DarkTurquoise' 'DarkOrange' 'Blue' 'Crimson'}); %jet(nsim);
figure
hax = gca;
hzoom = axes('position',[0.6445    0.1727    0.3205    0.3172]);
subplot(hax),  hold on
for isim=1:nsim
    % load results
    load(fullfile(local,outputpathref{isim},fmecadbfile)) % reload variable fmecadb
    nodes = fieldnames(fmecadb)';
    % build parenting tree
    for imol = 1:nmol % for each molecule
        if imol==1
            molnodes   = uncell(regexp(nodes,['^' mollist{imol} '.*'],'match'),[],[],true)'; % all nodes starting by X, Y
            nnodes = length(molnodes); % nodes
            CF =  cellfun(@(p) fmecadb.(p).CF,molnodes)'; % parents
            SML = cellfun(@(p) fmecadb.(p).SML,molnodes)'; % parents
        else
            molnodes = regexprep(molnodes,sprintf('^%s',mollist{imol-1}),mollist{imol});
            CF(:,imol) =  cellfun(@(p) fmecadb.(p).CF,molnodes)'; % parents
            SML(:,imol) = cellfun(@(p) fmecadb.(p).SML,molnodes)'; % parents
        end
        parents = cellfun(@(p) fmecadb.(p).parent,molnodes,'UniformOutput',false); % parents
        [~,~,~,childrens] = buildmarkov(molnodes,parents); % list column-wise child indices of each node
        adj = sparse(nnodes,nnodes); for i=1:nnodes, adj(i,childrens(childrens(:,i)>0,i))=1; end % adjency matrix
    end
    % graph
    s = severity(CF,SML);
    gplot(adj,s,'o-')
end
% update colors
hp=flipud(get(hax,'children')); for isim=1:nsim, set(hp(isim),'Color',col(isim,:),'linewidth',2); end % apply colors
% limit values
plot(xlim,[1 1]*100,'k:','linewidth',2); %plot(xlim,[1 1]*SML(1,2),'k:','linewidth',2)
plot([1 1]*100,ylim,'k:','linewidth',2); %plot([1 1]*SML(1,1),ylim,'k:','linewidth',2)
xlabel('severity (X)','fontsize',16); %xlabel('C_X (kg\cdotm^{-3})','fontsize',16)
ylabel('severity (Y)','fontsize',16); %ylabel('C_Y (kg\cdotm^{-3})','fontsize',16)
hl = legend(hp,{'design A0' 'design B0' 'design A1' 'design B1'},'location','NorthEastOutside','fontsize',14); set(hl,'box','off')
formatax(hax,'fontsize',14), axis tight
plotdata = gcfd;
% add zoom
subplot(hzoom), scfd(plotdata(1),'noaxes','nolegend'); axis([-0.2 5 -0.2 5])
formatax(hzoom,'fontsize',12)
% print
print_pdf(600,figname,figurefolder,'nocheck')
print_png(200,figname,figurefolder)

 
%% interpretation based on severity (for checking only)
% severity=@(CF,SML) 100*99./max(100*SML./CF-1,0); % Severity definition
% isim = 3;
% load(fullfile(local,outputpathdel{isim},fmecadbfile)) % reload variable fmecadb
% [res,resunique] = fmecasingle(fmecadb); % contribution of single step to CF
% %contribution of a single step as severity
% ress = cell2struct( cellfun( @(step,sml) severity(res.(step),sml),...
%                              fieldnames(res),...
%                              cellfun(@(step) fmecadb.(step).SML,fieldnames(res),'UniformOutput',false),...
%                              'UniformOutput',false ...
%                            ), ...
%                      fieldnames(res) );



%% Severity as a PARETO chart
severity=@(CF,SML) 100*99./max(100*SML./CF-1,0); % Severity definition
colinterp = @(s) interp1(linspace(0,1,64),jet(64),s/100,'linear',0); % Severity color (from blue to red and black beyond 100)
figname = [regexprep(fmecamainfile,'\..*$','') '_pareto']; % figure name
hfig=figure; formatfig(hfig,'figname',figname,'paperposition',[3.6154    2.0887   13.7533   25.5000 ]); % page format
hs = subplots([1 1],[ones(1,nsim) .5],.03,[ones(1,nsim-1)*.15 0],'alive',[1:4 6:9]); %subplots([.2 1],[.6 .4],[],[],'alive',3); 
hs = reshape(hs,nsim,2);
leg = {'A0' 'B0' 'A1' 'B1'};
for isim=1:nsim
    load(fullfile(local,outputpathdel{isim},fmecadbfile)) % reload variable fmecadb
    [res,resunique] = fmecasingle(fmecadb);
    nodes = fieldnames(res);
    ymax = 0;
    for imol = 1:nmol % for each molecule
        molnodes   = uncell(regexp(nodes,['^' mollist{imol} '.*'],'match'),[],[],true)'; % all nodes starting by X, Y
        nnodes = length(molnodes); % nodes
        subplot(hs(isim,imol)), hold on, i = 0;
        for step=molnodes % for all nodes stored in res
            i=i+1; % counter
            s   = severity(res.(step{1}),fmecadb.(step{1}).SML); % severity for current node
            col = colinterp(s);
            s(s==0)=NaN; % force NaN for severities=0
            for j=1:length(s), plot(i,s(j),'markerfacecolor',col(j,:),'markeredgecolor',col(j,:),'linestyle','none','markersize',8,'marker','o'), end
            text(i,0,[step{1} '  '],'fontsize',9,'rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle');
        end
        set(hs(isim,imol),'xlim',[0 nnodes+1])
        for r = [33 50 100]
            hr = refline(0,r); set(hr,'color',colinterp(r),'linewidth',2,'linestyle',':');
            if imol==2
                text(max(xlim),r,sprintf('Severity=%d',r),'HorizontalAlignment','right','VerticalAlignment','bottom','fontsize',8,'fontweight','bold','color',colinterp(r));
            end
        end
        xax = xlim; axis tight; xlim(xax)
        ymax = max(ymax,max(ylim));
        if imol==1, ylabel(sprintf('Severity %s',leg{isim}),'fontsize',14); end
    end
    set(hs(isim,:),'ylim',[0 ymax]), 
end
formatax(hs,'fontsize',12,'xticklabel',''); %% set(hs(:,2),'yticklabelmode','auto')
titles(hs,[],'x',.9,'y',.9,'suffix',')','fontsize',12)
%print
print_png(300,figname,figurefolder);
print_pdf(600,figname,figurefolder,'nocheck');
