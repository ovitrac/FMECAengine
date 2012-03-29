function rankout = paretochart(values,nodenames,valmax,defaultlayout,xlabeltxt,varargin)
% PARETOCHART: generate pareto chart (like in fmecaengine.m)
% SYNTAX: rankout = paretochart(values,nodesnames,valmax,defaultlayout,xlabeltxt,[property1,value1,property2,value2,...])
%        values: values to rank
%     nodenames: names correspond to values (label for axis y)
%        valmax: maximum value for color scale
% defaultlayout: true for default layout (main plot)
%     xlabeltxt: label for axis "x" default empty
%      titletxt: alternative text
%
%   Pair properties: Any valid axes property is accepted, in addition with
%       'xtxt': position for text (default = [], automatic)
%   'colormap': full colormap (default = jet(64))
%      Note use xscale': 'log' to display a log scale instead of a linear scale
% 
% Migration 2.0 - 13/01/12 - INRA\Olivier Vitrac - rev. 29/03/12

% Revision history
% 08/02/12 add xtxt,add axes properties
% 29/03/12 add colormap as property, accept [] as valax and defaultlayout

% default
default = struct('xtxt',[],'colormap',jet(64));
defaultax = struct('xscale','linear');
defaultcol = [0 0 0];

% arg check
if nargin<1, error('1 inputs are required'), end
if nargin<3, valmax = []; end
if nargin<4, defaultlayout = []; end
if nargin<5, xlabeltxt = ''; end
if isempty(valmax), valmax = max(values); end
if isempty(defaultlayout), defaultlayout = true; end

% if nargin<6, titletxt = ''; end
[options,optionsax] = argcheck(varargin,default);
optionsax = argcheck(optionsax,defaultax,'','keep');

% colors
ncol = size(options.colormap,1);
if size(options.colormap,2)~=3 || ncol<5, error('colormap must be a nx3 array with n>5'); end

nodenames = regexprep(nodenames,'_','\\_'); % protect '_'
[values,rank] = sort(values,'ascend'); %locvaluesmax = max(values);
col = interp1(linspace(0,1,ncol),options.colormap,values/valmax);
if any(isnan(col)), col(any(isnan(col),2),:) = 0; end
hold on
if isempty(options.xtxt)
    switch lower(optionsax.xscale)
        case 'linear', options.xtxt = min(0,0 - 0.01 * (max(values)-min(values)));
        case 'log',    options.xtxt = 0.9;
        otherwise, error('unknown xscale property ''%s'', only ''linear'' or ''log'' are accepted',optionsax.xscale)
    end
end
for inode = 1:length(nodenames)
    barh(inode,values(inode),'FaceColor',col(inode,:))
    labeltxt = nodenames(rank(inode));
    if defaultlayout
        text(0,inode,[labeltxt ' '],'fontsize',8,'HorizontalAlignment','right','VerticalAlignment','middle')
        text(values(inode),inode,sprintf(' %0.3g',values(inode)),'fontsize',10,'HorizontalAlignment','left','VerticalAlignment','middle')
    else
        if values(inode)>options.xtxt, coltxt = defaultcol;
        elseif sum(col(inode,:))/3<0.4, coltxt = rgb('Ivory');
        else coltxt = defaultcol;
        end
        text(options.xtxt,inode,labeltxt,'fontsize',6,'HorizontalAlignment','right','VerticalAlignment','middle','FontWeight','normal','color',coltxt)
    end
end
if defaultlayout, xlabel(xlabeltxt,'fontsize',16), end
% if isempty(titletxt), title(strrep(sprintf('%s[\\bf%s\\rm]',o.fmecamainfile,o.fmecasheetname),'_','\_'),'fontsize',12)
% else title(titletxt,'fontsize',10)
% end
axis tight, set(gca,'yticklabel',' ',optionsax)
if nargout, rankout = rank; end