function hsout = printtable(tr,th,styletr,styleth,stylesheet,varargin)
% PRINTTABLE draws table with predefined contents
% INPUTS
%               tr: n x m cell containing contents to put in table (n rows and m columns), data accepting: string, number, image... 
%               th: l x m cell containing header contents (l rows and m column)
%          styletr: style of contents
%          styleth: style of headers
%       stylesheet: structure defining styles
% VARARGIN
%               parent: handle to draw table (default=[])
%        paperposition: position of figure, (default = [])
%            papertype: type of paper (A3,A4...), (default = 'A4')
%           fontsizeth: fontsize for text at headers (default = 10)
%           fontsizetd: fontsize for text at content (default = 8)
%            linewidth: line width of table (default=1)
%       nrowmaxperpage: max number of rows for each page (default = 20)
%horizontalalignmentth: horizontal alignement for text at headers (default = center)
%horizontalalignmenttd: horizontal alignement for text at content (default=center)
%  verticalalignmentth: vertical alignement for text headers (default=middle)
%  verticalalignmenttd: vertical alignement for text content (default=middle)
%         fontweightth: font weight for text at headers (default=bold)
%         fontweighttd: font weight for text at content (default=normal)
%            numlayout: anonymous function coding for numeric layout (default = [], no layout)
%                       examples: @(x) formatsci(x,'eco') 
%                                 @formatsci
%        tablesubplotx: x dimension for subplots (draw general table)(default=[])
%        tablesubploty: y dimension for subplots (draw general table)(default=[])
%          imgsubplotx: x dimension for subplots for image (if required) (default=[.2 .8 .1])
%          imgsubploty: y dimension for subplots for image (if required) (default=[.05 .9 .05])
%             imgalive: subplot to be kept for image (default=5)
%              printon: save figure (default=false)
%         figurefolder: figure to save figure (default=cd)
%           resolution: resolution of figure to be saved (default=300)
%
% OUTPUTS
%                    h: array of handles
%
% Examples
%{
    % full figure
    for i = 1:3, tr{i,1} = 'a'; tr{i,2}= '1';  tr{i,3}= 'm'; end
    printtable(tr,[],[],[],[],'paperposition',[0.4000    0.1774   20.5391   29.0000])
    % embedded figure
    figure, hs = subplots(1,[1 1],0,0); subplot(hs(2)),printtable(tr,[],[],[],[],'parent',hs(2))
%}

% MS v. 2.1 - 07/10/2013 - INRA\Olivier Vitrac, LNE\Mai Nguyen - rev. 29/07/17 


% Revision history
% 22/10/13 : add print in eps format
% 23/10/13 : add hsout
% 13/03/15 : fix empty inputs and update accordingly help
% 11/03/17 : fix the number of rows to match the number of headers
% 29/07/17 : add numlayout

% default
default = struct('parent',[],'nrowmaxperpage',20,'horizontalalignmentth','center','verticalalignmentth','middle','horizontalalignmenttd','center',...
                 'verticalalignmenttd','middle','figname','page','linewidth',1,'fontweightth','bold','fontweighttd','normal','fontsizeth',10,'fontsizetd',8,...
                 'numlayout',[],...
                 'paperposition',[],'tablesubplotx',[],'tablesubploty',[],'imgsubplotx',[.2 .8 .1],'imgsubploty',[.05 .9 .05],'imgalive',5,'printon',false,'figurefolder',cd,'resolution',300);
styletr_default = 'td';
styleth_default = 'th';

% argcheck
o = argcheck(varargin,default,'nostructexpand');
stylesheet_default = struct(... use HTML nomenclature
                                   'th',struct('x',.5,'y',.5,'HorizontalAlignment',o.horizontalalignmentth,'VerticalAlignment',o.verticalalignmentth,'FontWeight',o.fontweightth,'FontSize',o.fontsizeth),...
                                   'td',struct('x',.5,'y',.5,'HorizontalAlignment',o.horizontalalignmenttd,'VerticalAlignment',o.verticalalignmenttd,'FontWeight',o.fontweighttd,'FontSize',o.fontsizetd));

if nargin < 1, error('1 argument (table content) is required'), end
if nargin < 2, th = []; end
if ~isempty(th), if length(th)~=size(tr,2),error('Please check size of headers and contents, number of columns must be the same'), end, end

if nargin < 3, styletr = []; end
if nargin < 4, styleth = []; end
if nargin < 5, stylesheet = []; end

if isempty(styletr), styletr = styletr_default; end
if isempty(styleth), styleth = styleth_default; end
if isempty(stylesheet), stylesheet = stylesheet_default; end


if ~iscell(styletr), styletr = {styletr}; end    
if ~iscell(styleth), styleth = {styleth}; end

styletr(end+1:end+max(0,size(tr,2)-length(styletr))) = styletr(end);
styleth(end+1:end+max(0,size(th,2)-length(styleth))) = styleth(end);
   
fillcell = @(texte,style) text(style.x,style.y,texte,rmfield(style,{'x','y'})); % function to fill cell

% printting page-wise
nheaders = size(th,1); nrows = size(tr,1); ncolumn = size(tr,2);
currentrow = 0; currentpage = 0; %npages = max(1,ceil(nrows/o.nrowmaxperpage));
if ~isempty(o.tablesubplotx), tablesubplotx = o.tablesubplotx; else tablesubplotx = ones(1,ncolumn); end
if ~isempty(o.tablesubploty), tablesubploty = o.tablesubploty; else tablesubploty = ones(1,min(nrows,o.nrowmaxperpage)+nheaders); end
while currentrow<nrows
    % new page
    currentpage = currentpage + 1;
    if isempty(o.parent) 
        if ~isempty(o.figname), figname = o.figname; else, figname = sprintf('page%d',currentpage); end
        formatfig(figure,'figname',figname,'paperposition',o.paperposition)
    end
    hs = subplots(tablesubplotx,tablesubploty,0,0); validrow = true(min(nrows,o.nrowmaxperpage+nheaders),ncolumn);
%     title(sprintf('page \\bf%d\\rm/\\bf%d',currentpage,npages),'fontsize',8)
    % headers
    if ~isempty(th)
        for iheader = 1:nheaders
            for icol = 1:ncolumn
                subplot(hs(iheader,icol)), fillcell(th{iheader,icol},stylesheet.(styleth{icol}))
            end
        end
    end
    % content
    for irow = 1:min(nrows,o.nrowmaxperpage)
        currentrow = currentrow + 1;
        for icol = 1:ncolumn
            if currentrow<=size(tr,1)
                 if iscellstr(tr{currentrow,icol}) || ischar(tr{currentrow,icol}), subplot(hs(irow+nheaders,icol)), fillcell(tr{currentrow,icol},stylesheet.(styletr{icol}))
                 elseif isnumeric(tr{currentrow,icol}) && ~ismatrix(tr{currentrow,icol}), him = subplots(o.imgsubplotx,o.imgsubploty,0,0,'position',hs(irow+nheaders,icol),'alive',o.imgalive,'strict'); % image, insert a small axes inside to incorporate margin
                        subplot(him),imshow(tr{currentrow,icol})
                 elseif isnumeric(tr{currentrow,icol}), subplot(hs(irow+nheaders,icol))
                     if isempty(o.numlayout)
                         fillcell(tr(currentrow,icol),stylesheet.(styletr{icol}))
                     else
                         fillcell(o.numlayout(tr{currentrow,icol}),stylesheet.(styletr{icol}))
                     end
                 end                
            else
                validrow(irow+nheaders,:) = false;
                delete(hs(irow+nheaders,:))
            end
        end
    end
    % line layout
    formatax(hs,'xlim',[0 1],'ylim',[0 1],'xticklabel',' ','yticklabel',' ','xtick',[0 1],'ytick',[0 1],'linewidth',o.linewidth)
%     set(hs(validrow),'xlim',[0 1],'ylim',[0 1],'xticklabel',' ','yticklabel',' ','xtick',[0 1],'ytick',[0 1]), drawnow
    % print page per page
    if o.printon
        print_pdf(o.resolution,get(gcf,'filename'),o.figurefolder,'nocheck')
        print_png(o.resolution,get(gcf,'filename'),o.figurefolder,'',0,0,0)
        print_eps(o.resolution,get(gcf,'filename'),o.figurefolder)
    end
    hsout = hs;
end

