% SCRIPT TO ASSESS MODAL CRITICALITY of bolino case with "real" values
% INRA\Olivier Vitrac - 23/06/11

% execution  with fmecaengine v0.41

% rev : 
%% Local path to initiate the simulation (to be executed to assign variables)
switch localname
    case 'WSLP-OLIVIER2'
        local = 'C:\Data\Olivier\INRA\Etudiants & visiteurs\Audrey Goujon\illustration_fmeca\';
    case 'WS-MOL4'
        local = 'C:\Data\Olivier\Audrey_Goujon\Matlab_FMECA\production';
    case 'mol10.agroparistech.fr'
        local = '/home/olivier/MaiNguyen/ANR_SFPD/WP2_FMECA/modal_criticality';
end
outputpathref = fullfile('modalcriticality_bolino_case_1','reference');
outputpath = fullfile('modalcriticality_bolino_case_1','omission');
fmecadbfile = 'fmecamodal.mat';

%% Reference simulation to introduce 2x2 modalities: [t t+]Storage x [D D+]OvenHeating
fmecadef = struct( ...
'local',local,...
'fmecamainfile','bolino_case_1.ods',...
'inputpath','inputs',...
'outputpath',outputpathref ...
);
[fmecadb,rawdata,dataref] = fmecaengine(fmecadef); % dataref (with inherited fields) is keep for modifications

%% Simulations with omission
data = dataref;
[data(1:length(data)).Foscale] = deal('[0 0]'); % add deletion "code" to all rows
fmecaengine(...
'local',local,...
'fmecamainfile',data,...
'inputpath','omission.ods',...
'outputpath',outputpath, ...
'fmecadbfile',fmecadbfile ...
);

%% Interpretation CF(path length = 4) - CF(path length = 3)
load(fullfile(local,outputpath,fmecadbfile)) % load fmecadb (the variable is created)
[res,resunique] = fmecasingle(fmecadb);

%% PARETO
severity=@(CF,SML) 100*99./max(100*SML./CF-1,0); % Severity definition
colinterp = @(s) interp1(linspace(0,1,64),jet(64),s/100,'linear',0); % Severity color (from blue to red and black beyond 100)
hfig=figure; formatfig(hfig,'figname',sprintf('PARETO_%s',regexprep(fmecadbfile,'.mat','')),'paperposition',[ 5.0467    7.2245   10.8907   15.2284]); % page format
hs = subplots([.2 1],[.6 .4],[],[],'alive',3); % layout 
subplot(hs), hold on, i = 0;
for step=fieldnames(res)' % for all nodes stored in res
    i=i+1; % counter
    s   = severity(res.(step{1}),fmecadb.(step{1}).SML); % severity for current node
    col = colinterp(s);
    s(s==0)=NaN; % force NaN for severities=0
    for j=1:length(s), plot(i,s(j),'markerfacecolor',col(j,:),'markeredgecolor',col(j,:),'linestyle','none','markersize',12,'marker','o'), end
    text(i,0,[step{1} '  '],'fontsize',16,'rotation',90,'HorizontalAlignment','right','VerticalAlignment','middle');
end
for r = [33 50 100]
    hr = refline(0,r); set(hr,'color',colinterp(r),'linewidth',2,'linestyle',':');
    text(max(xlim),r,sprintf('Severity=%d',r),'HorizontalAlignment','right','VerticalAlignment','bottom','fontsize',10,'fontweight','bold','color',colinterp(r));
end
formatax(hs,'fontsize',14,'xticklabel','','xlim',[0 i+1]), ylabel('Severity','fontsize',16);
print_png(300,get(gcf,'filename'),fullfile(local,outputpath));

%% RISK RANKING
steplist = fieldnames(resunique.value); nsteps = length(steplist);
colinterp = @(s) interp1(linspace(0,1,64),jet(64),s/100,'linear',0); % Severity color (from blue to red and black beyond 100)
hfig=figure; formatfig(hfig,'figname',sprintf('RISK_%s',regexprep(fmecadbfile,'.mat','')),'paperposition',[  5.0467    5.0887   10.8907   19.5000]); % page format
hscopy = subplots(1,ones(1,nsteps),0,0.05); set(hscopy,'visible','off');
hs = subplots(1,ones(1,nsteps),0,0.05,'position',gcf);
smax = max(cellfun(@(n) max(severity(resunique.value.(n),fmecadb.(n).SML)),steplist)); % maximal severity
smin = min(cellfun(@(n) max(severity(resunique.value.(n),fmecadb.(n).SML)),steplist)); % minimal severity
set(hscopy,'xlim',[0 1.1*smax])
for i=1:nsteps
    dCF = resunique.value.(steplist{i}); % all dCF values
    [dCFu,iu] = unique(dCF); % unique values
    su = severity(dCFu,fmecadb.(steplist{i}).SML); % convert to severity
    proba = resunique.proba.(steplist{i}); % probability assiciated to dCF
    probau = arrayfun(@(u) sum(proba(dCF==u)), dCFu); %idem with dCFu
    ymax = 1.2*max(probau);
    % regions
    subplot(hscopy(i))
    colormap(rgb({'MintCream' 'MistyRose'}));  [xp,yp] = meshgrid([33 100 max(110,1.1*smax)],[0 ymax]);
    hp=patch('Vertices',[xp(:) yp(:)],'Faces',[1 2 4 3 1;3 4 6 5 3],'facecolor','flat','edgecolor','none','FaceVertexCData',[1;2],'CDataMapping','direct');
    set(hp,'visible','on'), ylim([0 ymax])
    % stem
    subplot(hs(i)), hold on
    col = colinterp(su); % interpolated colors
    for j=1:length(su)
        plot(su(j)*[1 1],[0 probau(j)],'color','k','linestyle','--','linewidth',2.5)
        plot(su(j),probau(j),'marker','o','markersize',12,'markerfacecolor',col(j,:),'markeredgecolor',col(j,:))
        text(su(j),0,[resunique.id.(steplist{i}){iu(j)} ' '],'fontsize',8,'HorizontalAlignment','right','VerticalAlignment','Bottom','Rotation',270)
    end
    ylim([0 ymax])
    ylabel(steplist{i},'fontsize',14)
    % references values
    for r = [10 33 50 100]
        plot([r r],ylim,'color',colinterp(r),'linewidth',2,'linestyle',':');
        text(r,max(ylim),sprintf('%d',r),'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',10,'fontweight','bold','color',colinterp(r),'rotation',270);
    end
end
formatax(hs,'fontsize',14,'xticklabel','','xlim',[0 1.1*smax],'xscale','linear','color','none')
set(hs(end),'xticklabelmode','auto')
xlabel('Severity','fontsize',16);
print_png(300,get(gcf,'filename'),fullfile(local,outputpath));

%% RISK RANKING
steplist = fieldnames(resunique.value); nsteps = length(steplist);
colinterp = @(s) interp1(linspace(0,1,64),jet(64),s/100,'linear',0); % Severity color (from blue to red and black beyond 100)
hfig=figure; formatfig(hfig,'figname',sprintf('RISK_%s',regexprep(fmecadbfile,'.mat','')),'paperposition',[  5.0467    5.0887   10.8907   19.5000]); % page format
hscopy = subplots(1,ones(1,nsteps),0,0.05); set(hscopy,'visible','off');
hs = subplots(1,ones(1,nsteps),0,0.05,'position',gcf);
smax = max(cellfun(@(n) max(severity(resunique.value.(n),fmecadb.(n).SML)),steplist)); % maximal severity
smin = min(cellfun(@(n) max(severity(resunique.value.(n),fmecadb.(n).SML)),steplist)); % minimal severity
set(hscopy,'xlim',[0 1.1*smax])
for i=1:nsteps
    dCF = resunique.value.(steplist{i}); % all dCF values
    [dCFu,iu] = unique(dCF); % unique values
    su = severity(dCFu,fmecadb.(steplist{i}).SML); % convert to severity
    proba = resunique.proba.(steplist{i}); % probability assiciated to dCF
    probau = arrayfun(@(u) sum(proba(dCF==u)), dCFu); %idem with dCFu
    ymax = 1.2*max(probau);
    % regions
    subplot(hscopy(i))
    colormap(rgb({'MintCream' 'MistyRose'}));  [xp,yp] = meshgrid([33 100 max(110,1.1*smax)],[0 ymax]);
    hp=patch('Vertices',[xp(:) yp(:)],'Faces',[1 2 4 3 1;3 4 6 5 3],'facecolor','flat','edgecolor','none','FaceVertexCData',[1;2],'CDataMapping','direct');
    set(hp,'visible','on'), ylim([0 ymax])
    % stem
    subplot(hs(i)), hold on
    col = colinterp(su); % interpolated colors
    for j=1:length(su)
        plot(su(j)*[1 1],[0 probau(j)],'color','k','linestyle','--','linewidth',2.5)
        plot(su(j),probau(j),'marker','o','markersize',12,'markerfacecolor',col(j,:),'markeredgecolor',col(j,:))
        text(su(j),0,[resunique.id.(steplist{i}){iu(j)} ' '],'fontsize',8,'HorizontalAlignment','right','VerticalAlignment','Bottom','Rotation',270)
    end
    ylim([0 ymax])
    ylabel(steplist{i},'fontsize',14)
    % references values
    for r = [10 33 50 100]
        plot([r r],ylim,'color',colinterp(r),'linewidth',2,'linestyle',':');
        text(r,max(ylim),sprintf('%d',r),'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',10,'fontweight','bold','color',colinterp(r),'rotation',270);
    end
end
formatax(hs,'fontsize',14,'xticklabel','','xlim',[0 1.1*smax],'xscale','linear','color','none')
set(hs(end),'xticklabelmode','auto')
xlabel('Severity','fontsize',16);
print_png(300,get(gcf,'filename'),fullfile(local,outputpath));


%% Severity GRAPH restricted to resunique (all combined)
res_asseverity = resunique.value; for f=fieldnames(resunique.value)', res_asseverity.(f{1}) = severity(resunique.value.(f{1}),fmecadb.(f{1}).SML); end
[hg,hgtmp]=fmecagraph(fmecadb,res_asseverity,'min',0,'max',100); delete(hgtmp)
set(hg,'filename',sprintf('SEVERITYGRAPH2_%s',regexprep(fmecadbfile,'.mat','')))
graphfilename = fullfile(local,outputpath,[get(gcf,'filename') '.png']);
print_png(300,graphfilename); pngtruncateim(graphfilename,true,100);
print_pdf(300,regexprep(graphfilename,'\.png$','.pdf'),filesep,'nocheck')
