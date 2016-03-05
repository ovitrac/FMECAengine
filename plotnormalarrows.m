function hpout = plotnormalarrows(x,y,xi,varargin)
%PLOTNORMALARROWS plot arrows above or below a curve (x,y) at points xi (the curve must be already plotted at the proper scale)
% Syntax: hp = plotnormalarrows(x,y,xi ['property1',value1,'property2',value2,...])
% List of pair property/value (ARROW(Start,Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir))
%                                                               B
%                                                              /|\           ^
%      length:  length of the arrowhead in pixels             /|||\          |
%   baseangle:  base angle in degrees (ADE)                  //|||\\        L|
%    tipangle:  tip angle in degrees (ABC).                 ///|||\\\       e|
%       width:  width of the base in pixels                ////|||\\\\      n|
%        page:  Use hardcopy proportions                  /////|D|\\\\\     g|
%                                                        ////  |||  \\\\    t|
%   autoscale: 'on' (default) for unequal axes          ///    |||    \\\   h|
%       scale: magnification/reduction factor          //<----->||      \\   |
%              (default value=1), use negative        /   base |||        \  V
%              value to get anti-normal arrows        E    angle||<-------->C
%                                                              |||tipangle
%       color: string or 1x3 rgb value (default='')            |||
%              combine 'facecolor' and 'edgecolor'             |||
%     faceolor: as above (default='')                          |||
%    edgecolor: as above (default='')                          |||
%                                                              |||
%         text: string or cell array                           |||
%    textscale: text extension (default=1)                     |||
%verticalignment: middle                                       |||
%horzontalalignement: top                                      |||
%      fontsize: 10                                            |||
%    fontweight: normal                                        |||
%      fontname: arial                                         |||
%      rotation: 0                                             |||
%                                                              |||
%  See also: arrow                                          -->|A|<-- width
%
% Example:
%{
    t = linspace(-4*pi,4*pi,300)/10; y = 50*sin(10*t); plot(t,y,'k-','linewidth',2);
    plotnormalarrows(t,y,t(10:8:end-10),'length',5,'width',1,'tipangle',30,'color','b','scale',5);
    plotnormalarrows(t,y,t(10:8:end-10),'scale',-5,'length',5,'width',1,'tipangle',30,'color','r');
%}

% MS 2.1 - 26/02/2016 - INRA\Olivier Vitrac - rev. 29/02/2016

% Revision history
% 27/02/2016 RC
% 29/02/2016 add notunique

% default
default = struct(...
    'text','',...
    'textscale',1,...
    'length',[],...
    'baseangle',[],...
    'tipangle',[],...
    'width',[],...
    'page',[],...
    'scale',1,...
    'autoscale','on',...
    'notunique',false,...
    'yi',[]);
defaultcolor = struct(...
    'color','',...
    'facecolor','',...
    'edgecolor','' ...
    );
defaulttextproperties = struct(...
    'Horizontalalignment','center',...
    'Verticalalignment','middle',...
    'fontname','arial',...
    'fontsize',10,...
    'fontweight','normal',...
    'rotation',0,...
    'textcolor','k' ...
    );

% arg check
if nargin<3, error('hp = plotnormalarrows(x,y,xi [,''property1'',value1...)'), end
[o,remain] = argcheck(varargin,default);
addtext = ~isempty(o.text);
if ~iscell(o.text), o.text = {o.text}; end
if ~iscellstr(o.text), error('text must a be cell array of string or a string'), end 
[colors,remain] = argcheck(remain,defaultcolor);
textproperties = renfield(argcheck(remain,defaulttextproperties),'textcolor','color');
if ~isempty(colors.color)
    if isempty(colors.facecolor), colors.facecolor = rgb(colors.color); end
    if isempty(colors.edgecolor), colors.edgecolor = rgb(colors.color); end
end
colors = rmfield(colors,'color');
if isempty(colors.facecolor), colors = rmfield(colors,'facecolor'); end
if isempty(colors.edgecolor), colors = rmfield(colors,'edgecolor'); end
if o.notunique
    x = x(:); y = y(:);
else
    [x,ind] = unique(x(:)); y = y(ind); y = y(:);
end
xi = xi(:);
if ~isempty(o.yi)
    o.yi = o.yi(:);
    if length(o.yi) ~= length(xi), error('xi and yi must be the same length'), end
end
ni = length(xi);
drawnow
xax = xlim;
yax = ylim;

%% scaling
switch lower(o.autoscale)
    case 'on'
        dx = diff(xlim);
        dy = diff(ylim);
    case ''
    otherwise
        dx = 1;
        dy = 1;
end
x = (x-xax(1))/dx;
y = (y-yax(1))/dy;
xi = (xi-xax(1))/dx;

%% xi,yi
if isempty(o.yi)
    yi = interp1(x,y,xi,'pchip');
else
    yi = (o.yi-yax(1))/dy; 
end
dydx = ndf(x,y,1,[],'makeuniform',true);
if any(isinf(dydx))
    dyidxi = Inf(size(xi));
else
    dyidxi = interp1(x,dydx,xi,'pchip');
end
% figure, plot(xax(1)+x*dx,yax(1)+y*dy,'-',xax(1)+xi*dx,yax(1)+yi*dy,'.') %for development purposes
%% plot
% normals
% equation y -dyidxi*x  +dyidxi*xi-yi = 0
% normal = [1 -dyidxi]
fullscale = max(dx,dy);
hp = zeros(ni,1);
if addtext, o.text = argpad(o.text,ni); end
for i=1:ni
    translation = [-dyidxi(i) 1]; %[1 -dyidxi(i)];
    if isinf(translation(1))
        translation = o.scale/fullscale*[1 0];
    else
        translation = o.scale/fullscale*translation/norm(translation);
    end
    xa = xax(1) + (xi(i)+[0;translation(1)])*dx;
    ya = yax(1) + (yi(i)+[0;translation(2)])*dy;
    %line(xa,ya)
    hp(i)=arrow([xa(1) ya(1)],[xa(2) ya(2)],o.length,o.baseangle,o.tipangle,o.width,o.page);
    if addtext
        xt = xax(1) + (xi(i)+translation(1)*o.textscale)*dx;
        yt = yax(1) + (yi(i)+translation(2)*o.textscale)*dy;
        text(xt,yt,o.text(i),textproperties)
    end
end

%% apply colors
if ~isempty(fieldnames(colors)), setprop(hp,colors), end

%% output
if nargout, hpout = hp; end