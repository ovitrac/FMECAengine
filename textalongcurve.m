function htout=textalongcurve(x,y,xt,txt,varargin)
%TEXTALONGCURVE display text tangentially along curves (a callback fix angles for best printing quality)
% LONG SYNTAX
%   ht=textalongcurve(x,y,xt,txt [,'property',value,...])
%   x,y : mx1 array setting the coordinates of the curve (one single curve)
%    xt : nx1 array setting the x-values where the text must be printed (y-values are interpolated)
%   txt : 1x1 or nx1 cell array containing the text ti display
%   property/value: any valid text pair properties or properties listed below
%       xscale: 'linear' or 'log'
%       yscale: 'linear' or 'log'
%    operation: anonymous function to define operation on angles (default=@(x) x), use @(x) 90+x to take normals instead of tangents
% SHORT SYNTAX
%   ht=textalongcurve(hp,xt,txt [,'property',value,...])
%   ht=textalongcurve(hax,hp,xt,txt [,'property',value,...])
%    hp : handle of a line such that x=get(hp,'Xdata'), y=get(hp,'Ydata'), hax=get(hp,'Parent')
%   hax : axes handle where to plot
%
%   NOTES: additional private properites are available in the code but not useful for most uses
%          xdir and ydir are not accounted in the present version
%          automatic xt assignment can be set using xt=[n n], where n is the number of text (size must be compatible)
%          angles are optimized for printing not for screens
%
%  %Example 1 - LONG SYNTAX (simple plot)
%       x = linspace(0,1.5,1e3)*pi; y=sin(x);
%       xt = linspace(0.25,1.25,5)*pi;
%       figure, plot(x,y), textalongcurve(x,y,xt,arrayfun(@(i) sprintf('test%d',i),1:length(xt),'UniformOutput',false))
%  %Example 1 - SIMPLE SYNTAX (simple plot)
%       x = linspace(0,1.5,1e3)*pi; y=sin(x);
%       xt = linspace(0.25,1.25,5)*pi;
%       figure, hp=plot(x,y); textalongcurve(hp,xt,arrayfun(@(i) sprintf('test%d',i),1:length(xt),'UniformOutput',false))
%
%  %Example 1bis - LONG SYNTAX (to highlight the capabilities of the callback)
%       figure, plot(x,y), textalongcurve(x,y,xt,'------------')
%  %Example 1bis - SHORT SYNTAX (to highlight the capabilities of the callback)
%       figure, textalongcurve(plot(x,y),xt,'------------')
%
%  %Example 1ter - SHORT SYNTAX (using normals)
%       figure, textalongcurve(plot(x,y),xt,'\rightarrow','operation',@(x)90+x,'horizontalalignment','left')   
%
%  %Example 2 - LONG SYNTAX (semi-log)
%       x = linspace(0,10,1e3); y=100*(1-exp(-.1*x));
%       figure, plot(x,y), set(gca,'yscale','log'), textalongcurve(x,y,logspace(-1,0.9,25),'test','yscale','log')
%  %Example 2 - SHORT SYNTAX (semi-log)
%       x = linspace(0,10,1e3); y=100*(1-exp(-.1*x));
%       figure, hp=plot(x,y); set(gca,'yscale','log'), textalongcurve(hp,logspace(-1,0.9,25),'test')
%
%  %Example 3 - LONG SYNTAX
%       x = logspace(1,3,1e3); y = x.^2;
%       figure, plot(x,y), set(gca,'xscale','log','yscale','log'), textalongcurve(x,y,[7 7],'tst','fontsize',8,'xscale','log','yscale','log')
%  %Example 3 - SHORT SYNTAX
%       x = logspace(1,3,1e3); y = x.^2;
%       figure, hp=plot(x,y); set(gca,'xscale','log','yscale','log'), textalongcurve(hp,[7 7],'tst','fontsize',8)

% MS 2.1 - 11/11/2012 - INRA\Olivier Vitrac - rev. 12/11/12

% Revision history
%12/11/12 add operation

% CALLBACK: special behavior for updating angles
if ~nargin, x=gcf; end
if numel(x)==1 && ishandle(x) && strcmpi(get(x,'Type'),'figure')
    fig = x;
    for child = get(fig,'children')'
        if strcmpi(get(child,'Type'),'axes')
            scale = scaleplot(fig,child);
            for htxt = findobj(child,'Tag','textalongcurve')'
                obj = get(htxt,'UserData');
                set(htxt,'rotation',obj.operation(atan(scale*obj.slope)*180/pi),'horizontalalignment',obj.horizontalalignment,'verticalalignment',obj.verticalalignment);
            end
        end
    end
    return
end

% SYNTAX with handles (short syntax)
if ishandle(x)
    if strcmpi(get(x,'Type'),'axes')
        if nargin<4, error('4 arguments are required'), end
        hax = x;
        if strcmpi(get(y,'Type'),'line')
            x = get(y,'Xdata');
            y = get(y,'Ydata');
        else
            error('the second argument must be the handle of a line object')
        end
    elseif strcmpi(get(x,'Type'),'line')
        if nargin<3, error('3 arguments are required'), end
        if nargin<4, tmp = []; else tmp = txt; end
        hax = get(x,'parent');
        txt = xt;
        xt = y;
        y = get(x,'Ydata');
        x = get(x,'Xdata');
        if ~isempty(tmp), varargin = [{tmp} varargin]; end
    else
        error('the first argument must the handle of a axes or line object');
    end
else
    if nargin<4, error('4 arguments are required: ht=textalongcurve(x,y,xt,txt [,property,value,...])'), end
    hax = gca;
end

% default
default = struct(...
    'method','cubic',...
    'ninterp',1e4,...
    'xscale',get(hax,'xscale'),...
    'yscale',get(hax,'yscale'),...
    'horizontalalignment','center',...
    'verticalalignment','middle',...
    'scale',[],...
    'nmin',11,...
    'operation',@(x) x ... added 12/11/12
    );

% arg check (common to short and long syntax)
[options,others] = argcheck(varargin,default,'','keep');
xax=get(hax,'xlim'); yax=get(hax,'ylim');
if strcmpi(options.xscale,'log'), x = log10(x(:)); xt=log10(xt(:)); xax=log10(xax); else x=x(:); xt=xt(:); end
if strcmpi(options.yscale,'log'), y = log10(y(:)); yax=log10(yax); else y=y(:); end
m = size(x,1); [xmin,xmax,ymin,ymax] = deal(xax(1),xax(2),yax(1),yax(2));
if size(y,1)~=m, error('x and y must be of the size'), end
ind = find(~isinf(x) & ~isnan(x) & ~isinf(y) & ~isnan(y));
if length(ind)<options.nmin, error('insufficient number of valid data (%d points are at least required)',options.nmin), end
if length(ind)<m, x = x(ind); y=y(ind); end
xmin = max(xmin,min(x)); xmax = min(xmax,max(x)); ymin = max(ymin,min(y)); ymax = min(ymax,max(y));
xt = xt(:); n = size(xt,1);
if n==2 && (xt(1)==xt(2))
    if strcmpi(options.xscale,'log'), n = 10.^xt(1); else n=xt(1); end
    xt = linspace(xmin+0.1*(xmax-xmin),xmax-0.1*(xmax-xmin),n);
end
if ~iscell(txt), txt = {txt}; end
if ~iscellstr(txt), error('txt must be a cell array'), end
if length(txt)==1, txt = repmat(txt,n,1); end
if length(txt)~=n, error('xt must be of the same size as xt'), end
yt = interp1(x,y,xt,options.method);
if isempty(options.scale)
    scale = scaleplot(gcf,gca);
end
    

% interpolation
xi = linspace(0,1,options.ninterp)';
yi = interp1((x-xmin)/(xmax-xmin),(y-ymin)/(ymax-ymin),xi,options.method);
slopet = interp1(xi,ndf(xi,yi),(xt-xmin)/(xmax-xmin),options.method);
alphat = atan(scale*slopet)*180/pi;
if strcmpi(options.xscale,'log'), xt =10.^xt; end
if strcmpi(options.yscale,'log'), yt =10.^yt; end

% text
ht = zeros(n,1);
% hold on, plot(xt,yt,'r+','markersize',12) % plot for debugging only
for i=1:n
    ht(i) = text(xt(i),yt(i),txt{i},'rotation',alphat(i),...
        'horizontalalignment',options.horizontalalignment,'verticalalignment',options.verticalalignment,others{:});
    set(ht(i),'Tag','textalongcurve',...
              'UserData',struct('slope',slopet(i),'operation',options.operation,...
                                'horizontalalignment',options.horizontalalignment,'verticalalignment',options.verticalalignment))
end

% output
set(gcf,'ResizeFcn',@textalongcurve) % callback
if nargout, htout = ht; end


%% function
function scale = scaleplot(fig,ax)
paperposition = get(fig,'paperposition'); position = get(ax,'position');
scale = paperposition(4)/paperposition(3)*position(4)/position(3);
