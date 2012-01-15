function rankout = paretochart(values,nodenames,valmax,defaultlayout,xlabeltxt)
% PARETOCHART: generate pareto chart (like in fmecaengine.m)
% SYNTAX: rankout = paretochart(values,nodesnames,[property1,value1,property2,value2,...])
%        values: values to rank
%     nodenames: names correspond to values (label for axis y)
%        valmax: maximum value for color scale
% defaultlayout: true for default layout (main plot)
%     xlabeltxt: label for axis "x" default "C_F"
%      titletxt: alternative text
% 
% Migration 2.0 - 13/01/12 - INRA\Olivier Vitrac - rev. 

% Revision history 


% arg check
if nargin<1, error('1 inputs are required'), end
if nargin<3, valmax = max(values); end
if nargin<4, defaultlayout = true; end
if nargin<5, xlabeltxt = ''; end
% if nargin<6, titletxt = ''; end
ncol = 64;

nodenames = regexprep(nodenames,'_','\\_'); % protect '_'
[values,rank] = sort(values,'ascend'); locvaluesmax = max(values);
col = interp1(linspace(0,1,ncol),jet(ncol),values/valmax);
if any(isnan(col)), col(any(isnan(col),2),:) = 0; end
hold on
for inode = 1:length(nodenames)
    barh(inode,values(inode),'FaceColor',col(inode,:))
    labeltxt = nodenames(rank(inode));
    if defaultlayout
        text(0,inode,[labeltxt ' '],'fontsize',10,'HorizontalAlignment','right','VerticalAlignment','middle')
        text(values(inode),inode,sprintf(' %0.3g',values(inode)),'fontsize',10,'HorizontalAlignment','left','VerticalAlignment','middle')
    else
        if sum(col(inode,:))/3<0.4 && values(inode)>.2*locvaluesmax, coltxt = rgb('Ivory'); else coltxt = [0 0 0]; end
        text(0,inode,labeltxt,'fontsize',7,'HorizontalAlignment','left','VerticalAlignment','middle','color',coltxt,'FontWeight','normal')
    end
end
if defaultlayout, xlabel(xlabeltxt,'fontsize',16), end
% if isempty(titletxt), title(strrep(sprintf('%s[\\bf%s\\rm]',o.fmecamainfile,o.fmecasheetname),'_','\_'),'fontsize',12)
% else title(titletxt,'fontsize',10)
% end
set(gca,'yticklabel',' '), axis tight
if nargout, rankout = rank; end