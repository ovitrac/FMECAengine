function [out,hsout,hfigout] = key2keybar(val,ttles,details,names,varargin)
% plot output of key2key request as pareto chart (deformulatio task of SFPD project)
% SYNTAX
%       key2keybar(val [,details/title,names,'property1',value1,'property2',value2,...)
%       [out,hsout,hfigout] = key2keybar(...)
%
% INPUTS:
%     val: output of key2key (assuming a nx1 cell array)  
% details: structure output of key2key or title
%   names: nx1 cell array giving as alternative names (depicted as a table below the pareto chart) 
%
% List of Properties/Values
%                'maxbar': maximum number of plotted bars (default = Inf)
%                'maxrow': maximum number of rows in the table below (default = Inf)
%           'rowheightcm': row height in cm (default=0.8)
%            'barwidthcm': bar width in cm (default=0.8)
%           'heightmincm': mheightmincminimum figure height in cm (default = 12)
%            'widthmincm': minimum figure width in cm (default = 15)
%               'figname': name of figure  (default='pareto')
%         'paperposition': paper position for figure (default=[])
%      'paperorientation': paperorientation (default='portrait')
%             'papertype': size, type of paper (default='A4)
%          'subplotspace': space between plot and table (default=0.1)
%             'rotationx': rotation angle for x ticklabel (default=90)
%           'fontsizeleg': fontsize for legend text in table (default=10)
%                 'title': text for title (default='')
%         'fontsizetitle': fontsize for title texte (default=16)
%             'fontsizey': fontsize for y label (default=16)
%  'horizontalalignmentx': horizontal alignment for x label (default='right')
%           'columntable': ratio of column in of the table below
%                  'cmap': colormap for bar plot (default=jet(64))
%            'resolution': resolution for saving figure (default=300)
%
% migration - 07/01/15 -  INRA\Olivier Vitrac, Mai Nguyen - rev. 18/8/2015
% history
% 19/1/2015: add axes for details/titles
% 18/8/2015: fix _ representation for details/titles

% default
papertype = struct('A4',[21 29.7],'A3',[29.7 42]);
default = struct('maxbar',Inf,'maxrow',Inf,'rowheightcm',0.8,'barwidthcm',0.8,'heightmincm',12,'widthmincm',15,'figname','pareto','paperposition',[],...
                 'paperorientation','portrait','papertype','A4','subplotspace',0.1,'fontsizex',12,'fontsizey',16,'fontsizeleg',10,'rotationx',90,...
                 'widthwraptext',40,'horizontalalignmentx','right','resolution',300,'columntable',[.3 .7],'title','','fontsizetitle',14,'tablealignment','center');
keyword = 'probab';

% argcheck
if nargin < 1, error('One argument is required'), end
if nargin < 2, ttles = {}; end
if nargin < 3, details = {}; end
if nargin < 4, names = {}; end
if ~iscellstr(val), error('val must be a cell array of strings'), end
if ~iscellstr(ttles), error('Titles must be a cell array of string'), end
if ~isstruct(details) && ~iscellstr(details), error('details must be a structure or cell array of string'), end
if isstruct(details), end    
if isempty(names), names = {}; end
if ~iscellstr(names),  error('names must be a cell array of strings'), end
val = val(:); n = length(val); if n==0; disp('Nothing to plot'), return, end
names = names(:); if ~isempty(names) && length(names)~=n, error('names must be the same size as val'), end
o = argcheck(varargin,default,keyword,'nostructexpand');

% remove empty entry
valid = arrayfun(@(m) ~isempty(val{(m)}),1:n);
val = val(valid);
if ~isempty(names), names = names(valid); end
% counts
[valu,ival,jval] = unique(val); 
counts = histc(jval,1:length(valu));
[countssort,isort] = sort(counts,'descend');
probabsort = countssort/sum(countssort);

% page layout
nbar = min(length(valu),o.maxbar);
if isempty(names), nrow = 0; else nrow = min(length(valu),o.maxrow); end
figwidth = max(nbar*o.barwidthcm,o.widthmincm);
figheight = max(nrow*o.rowheightcm,o.heightmincm);

papersize = papertype.(o.papertype);
if strcmpi(o.paperorientation,'landscape'), papersize = fliplr(papersize); end  
figwidth = min(papersize(1),figwidth);
if nrow > 0, figcoord = max(0,[(papersize(1)-figwidth)/2 (papersize(2)-figheight*2)/2]);
    if isempty(o.paperposition), o.paperposition = [figcoord figwidth  min(papersize(2)/2,figheight)*2]; end
else figcoord = max(0,[(papersize(1)-figwidth)/2 (papersize(2)-figheight)/2]);
    if isempty(o.paperposition), o.paperposition = [figcoord figwidth min(papersize(2)/2,figheight)]; end
end
hfig = figure; formatfig(hfig,'figname',o.figname,'paperposition',o.paperposition,'paperorientation',o.paperorientation,'papertype',o.papertype)
if o.probab, ytext = max(probabsort)/100; else ytext = max(counts)/100; end
% plot
if nrow>0
    hs = [subplots(1,[0.1 0.05*length(details) 1 1],0,[0.02 .02 o.subplotspace]);subplots(1,[0.1 0.05*length(details) 1 1],0,[0.02 .02 o.subplotspace],'alive',3,'position',gcf,'color','none','xAxisLocation','top','yAxisLocation','right')];
else
    hs = [subplots(1,[0.1 0.05*length(details) .9 .1],0,0,'alive',1:3); subplots(1,[0.1 0.05*length(details) .9 .1],0,0,'alive',3,'position',gcf,'color','none','xAxisLocation','top','yAxisLocation','right')];
end

subplot(hs(1)), text(0.5,0.5,regexprep(ttles,'_','\\_'),'fontsize',o.fontsizetitle,'horizontalalignment','center','verticalalignment','middle')
subplot(hs(2)), text(0.5,0.5,regexprep(details,'_','\\_'),'fontsize',8,'horizontalalignment','center','verticalalignment','middle')
formatax(hs(1:2),'xticklabel','','yticklabel','','visible','off')

subplot(hs(3)), hold on
for i = 1:length(counts)
    if o.probab, bar(i,probabsort(i),'facecolor',interp1(linspace(0,max(probabsort),64),jet(64),probabsort(i)))
    else bar(i,countssort(i),'facecolor',interp1(linspace(0,max(countssort),64),jet(64),countssort(i)))
    end
    text(i,-ytext,regexprep(valu(isort(i)),'_','\\_'),'rotation',o.rotationx,'horizontalalignment',o.horizontalalignmentx,'fontsize',o.fontsizex)
end
if o.probab, ylabel('Probability','fontsize',o.fontsizey), else ylabel('Counts','fontsize',o.fontsizey), end

subplot(hs(end)), hold on
if o.probab, plot(1:length(counts),cumsum(probabsort),'linewidth',1.5), ylabel('Cumulated probability','fontsize',o.fontsizey)
else plot(1:length(counts),cumsum(countssort),'linewidth',1.5), ylabel('Cumulated counts','fontsize',o.fontsizey)
end
formatax(hs([3 end]),'xlim',[0 length(valu)+1],'xticklabel','','fontsize',o.fontsizex,'box','off')

% table of leg
if nrow > 0
    namesu = names(ival); 
    txtintable = cell(length(valu),2);
    tabheight = zeros(1,length(valu));
    for i = 1:length(valu)
        txtintable{i,1} = regexprep(valu(isort(i)),'_','\\_') ; %valu{isort(i)};
        txtintable{i,2} = textwrap(namesu(isort(i)),o.widthwraptext);
        tabheight(i) = length(txtintable{i,2});
    end
    stylesheet = struct('td',struct('x',.5,'y',.5,'HorizontalAlignment',o.tablealignment,'VerticalAlignment','middle','FontWeight','normal','FontSize',o.fontsizeleg));
    styletr = 'td';
    subplot(hs(4)), printtable(txtintable,[],styletr,[],stylesheet,'nrowmaxperpage',length(valu),'linewidth',1,'parent',hs(4),'tablesubplotx',o.columntable,'tablesubploty',tabheight);
end
% output
if nargout
    if nrow==0
        out = struct('val',{valu(isort)},'names',{{}},'counts',countssort);
    else
        out = struct('val',{valu(isort)},'names',{namesu(isort)},'counts',countssort);
    end
    if nargout>1, hsout = hs; end
    if nargout>2, hfigout = hfig; end
end