function [hlout,htout] = lineannotation(xy,txt,varargin)
%LINEANNOTATION add an annotation line to current plot
%
%   Syntax: [hl,ht] = lineannotation(xy,txt [,'property',value,...])
%       xy: nx2 array as [x y] points
%      txt: nx1 cell array
%
%   List of pair property/values
%   ============================
%    shape: mx2 array defining a line as [0 0;[dx dy]] (default = [0 0; 1 1 ; 1 0)])
%           dx and dy are dimensionless displacements to create the line [x(i)+dx y(i)+dy] for all i=1..n
%
%    scale: scalar ranging between 0 and 1 (default = 0.05) setting the size of annotation lines relatively to axes
%
%   orientationx, orientationy are used to change the orientation of dx and dy (default value = 1)
%
% Text properties (see text for details, only default values are provided)
%       dxtxt: 0.0100 (text margin along x)
%       dytxt: 0.0100 (text margin along y)
%    fontname: 'arial'
%    fontsize: 8
%  fontweight: 'normal'
%   fontangle: 'normal'
%    rotation: 0
% verticalalignment: 'middle'
%   textcolor: 'k'
%
% Line properties (see line for details, only default values are provided)
%   linestyle: '-'
%   linewidth: 0.5000
%       color: 'k'
%
% Arrow properties (to be used along with keyword 'arrow', see arrow for details)
%      length: 10
%   baseangle: 30
%    tipangle: 40
%       width: 0.5
%
% example:
%   x=linspace(0,6*pi,1e3)'; y=sin(x);
%   x0=(pi/2:pi:x(end))'; y0 = sin(x0);
%   plot(x,y,'b-',x0,y0,'rx')
%   txt = cellfun(@(i) sprintf('peak %d',i),num2cell(1:length(x0),1),'UniformOutput',false);
%   lineannotation([x0 y0],txt,'orientationx',[1 -1 1 -1 1 -1],'orientationy',[1 -1 1 -1 1 -1]*2,'arrow')

% INRA\MS 2.1 - 21/05/2013 - Olivier Vitrac - rev. 25/05/2013

% Revision history
% 23/05/2013 release candidate
% 25/05/2013 enable text on several lines when n=1

% default
keyword = 'arrow';
default = struct(...
    'shape',[1 1; 1 0],... 0 0 for the starting point
    'scale',0.05,...
    'orientationx',1,...
    'orientationy',1,...
    'dxtxt',0.03,...
    'dytxt',0,...
    'linestyle','-',...
    'linewidth',0.5,...
    'color','k',...
    'fontsize',8,...
    'fontweight','normal',...
    'fontangle','normal',...
    'fontname','arial',...
    'textcolor','k',...
    'rotation',0,...
    'verticalalignment','middle',...
    'length',10,...
    'baseangle',30,...
    'tipangle',40,...
    'width',0.5 ...
    );

% arg check
if nargin<2, error('two arguments are required'), end
[n,nc] = size(xy);
if nc>2 && n==2, xy = xy'; [n,nc] = size(xy); end
if (nc~=2) || ~isnumeric(xy), error('xy must be a nx2 array'), end
if ischar(txt), txt = {txt}; end
if length(txt)==1, txt = txt(ones(n,1)); end
if length(txt)~=n
    if n==1, txt={txt}; else error('txt and xy must be of the same size'); end
end
o = argcheck(varargin,default,keyword);
[m,mc] = size(o.shape);
if m==2 && mc>2, o.shape=o.shape'; [m,mc] = size(o.shape); end
if (mc~=2) , error('shape must be a mx2 array'), end
o.orientationx = o.orientationx(:); nox = length(o.orientationx);
o.orientationy = o.orientationy(:); noy = length(o.orientationy);
o.orientationx = [o.orientationx;o.orientationx(nox*ones(max(0,n-nox),1))];
o.orientationy = [o.orientationy;o.orientationy(noy*ones(max(0,n-noy),1))];

% axes properties
drawnow
xgca = get(gca,'xlim');
ygca = get(gca,'ylim');
xscale = diff(xgca) * o.scale;
yscale = diff(ygca) * o.scale;

% build line array
[lx,ly] = deal(zeros(m+1,n));
lx(1,:) = xy(:,1)'; ly(1,:) = xy(:,2)';
for i=1:n % for all annotations
    for j=1:m % for all segments
        lx(j+1,i) = lx(j,i) + o.orientationx(i) .* o.shape(j,1)*xscale;
        ly(j+1,i) = ly(j,i) + o.orientationy(i) .* o.shape(j,2)*yscale;
    end
end

% arrow tip if needed
if m>1 && o.arrow
    for i=1:n
        ha = arrow([lx(2,i) ly(2,i)],[lx(1,i) ly(1,i)],o.length,o.baseangle,o.tipangle,o.width);
        set(ha,'facecolor',o.color,'edgecolor',o.color) % force the same color as lines
    end
    lx = lx(2:end,:); ly = ly(2:end,:); % remove the first segment
end

% plot lines
hl = line(lx,ly,'linestyle',o.linestyle,'linewidth',o.linewidth,'color',o.color);

% text
ht = zeros(n,1);
for i=1:n % for all annotations
    if o.orientationx(i)>0, horizontalalignment = 'left'; else horizontalalignment = 'right'; end
    ht(i) = text(...
         lx(end,i) +  o.orientationx(i) .* o.dxtxt*xscale,...
         ly(end,i) +  o.orientationy(i) .* o.dytxt*yscale,...,...
         txt{i},'horizontalalignment',horizontalalignment,...
                           'fontsize',o.fontsize,...
                         'fontweight',o.fontweight,...
                           'fontname',o.fontname,...
                           'rotation',o.rotation,...
                  'verticalalignment',o.verticalalignment,...
                              'color',o.textcolor ...
         );
end

% optional outputs
if nargout, hlout = hl; end
if nargout>1, htout = ht; end