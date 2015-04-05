function hfigout = formatfig(hfig,varargin)
%FORMATFIG format figures for publications (figname is used to set automatically the name and the filename)
%   syntax: formatfig
%           formatfig(hfig,'figname',figname,'property1',value1,'property2',value2,'property3',value3,..)
%           hfig = formatfig(...)
%      hs: axes handles (default = all handles in the current figure)
%      property/value any valid axes properties
%      property 'figname' codes for the name of the file (' ' and ':' are replaced respectively by '_' and '-')

% MS 2.1 - 28/12/09  - INRA\Olivier Vitrac - rev. 12/03/15

% Revision history
% 22/04/11 force case sensitive values (required for filename on Linux)
% 04/09/11 replace ifig by hfig(ifig) for figure numbering
% 08/01/15 new color scheme for defaultAxesColorOrder based R2014b http://fr.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
% 28/02/15 add Tag
% 12/03/15 add hfigout

arg_default = struct(...
'figname'  ,sprintf('fig_%s',datestr(now)),...
'Tag','',...
'PaperUnits','Centimeters',...
'PaperType','A4',...
'PaperOrientation','Portrait',...
'PaperPositionMode','manual',...
'Units','pixels',...
'defaultAxesColorOrder',[
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
] ...
);


% arg check
if nargin<1, hfig = []; end
if isempty(hfig)
    children = get(0,'children');
    hfig = children(strcmp(get(children,'type'),'figure'));
    hfig = hfig(end:-1:1);
end
if isempty(hfig), hfig = gcf; end
[arg,remain] = argcheck(varargin,arg_default,[],[],'case');
arg.figname(arg.figname==' ')='_';
arg.figname(arg.figname==':')='-';

% set attributes to figure(s)
set(hfig,rmfield(arg,'figname'))
if ~isempty(remain)
    try set(hfig,remain{:})
    catch, error('invalid additional (figure) arguments')
    end
end
nfig = length(hfig);

% change the name in the titlebar
if ~isempty(arg.figname)
    for ifig=1:nfig
        % before 4/9/11: set(hfig(ifig),'Name',sprintf('%0.2d: %s',ifig,arg.figname),'FileName',arg.figname,'NumberTitle','off')
        set(hfig(ifig),'Name',sprintf('%0.2d: %s',hfig(ifig),arg.figname),'FileName',arg.figname,'NumberTitle','off')
    end
end

% output
if nargout, hfigout = hfig; end