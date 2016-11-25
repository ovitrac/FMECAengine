function [bout,statout]=autoreglin(x,y,minx,alpha,ploton)
%AUTOREGLIN automatic search of best local linear regression lines (sliding bi-cubic weighting kernel)
%   syntax: [b,stat]=autoreglin(x,y[,minx,alpha,ploton])
%       x: mx1 vector
%       y: mx1 vetor such that y = f(x)
%    minx: minimum number of points to be included (default = 10)
%   alpha: probability (default = 0.05)
%  ploton: plot flax (default = true)

% IGA 1.0 - 21/11/08 - INRA\Olivier Vitrac - rev. 04/06/2016 

% Revision history
% 02/12/08 add minx
% 22/12/08 fix column vectors
% 04/08/10 use Xdata and Ydata as inputs
% 04/06/16 add residues and iok to stats (to facilitate bootstraping)

% Default
alpha_default = 0.05;
ngrid = 50;
minx_default = 10;   % minimun numbr of data to be included
nspline = 2; % criteria takes intoccount compare with data nspline times outside the initial window
maxiter = 10;
reltol = 0.001;
maxplot = 50;
ploton_default = true;
autoattrib = false;

% arg check
if ~nargin
    childs = get(gca,'children');
    if any(childs), x = get(childs(1),'Xdata'); y = get(childs(1),'Ydata'); autoattrib=true; end
end
if nargin<2 && ~autoattrib, error('syntax: [b,stat]=autoreglin(x,y[,x0,alpha])'), end
if nargin<3, minx = []; end
if nargin<4, alpha = []; end
if nargin<5, ploton = []; end
x = x(:); y = y(:);
i0ok =find(~isnan(y)); x = x(i0ok); y = y(i0ok);
xm = mean(x); ym=mean(y); sx = std(x); sy = std(y);
xn = (x-xm)/sx; yn = (y-ym)/sy; dxn = abs(diff(xn([1 end]))); dxnm = abs(mean(diff(xn)));
if isempty(minx), minx = minx_default; end
if isempty(alpha), alpha = alpha_default; end
if isempty(ploton), ploton = ploton_default; end

% Definitions
o = ones(size(x));
vy = var(y);


% Optimization
it = 0; relerr = Inf; 
xn_explore  = xn([1 end]);
dxn_explore = [max(minx*dxnm,dxn/ngrid) dxn];
F10last = [];
i10last = [];
Fminlast = +Inf;
while (relerr>reltol) && (it<maxiter)
    it = it +1 ;
    dispf('AUTOREGLIN step %d of %d (max)...',it,maxiter), t = clock;
    [X0,H] = meshgrid(linspace(xn_explore(1),xn_explore(2),ngrid),linspace(dxn_explore(1),dxn_explore(2),ngrid));
    F = arrayfun(@F_,X0,H);
    F10 = prctile(F(:),10);
    i10 = find(F<=F10);
    xn_explore  = [min(X0(i10)) max(X0(i10))];
    dxn_explore  = [min(H(i10)) max(H(i10))];
    Fmin = min(F(:));
    if Fmin<Fminlast
        relerr = 1-Fmin/Fminlast;
        Fminlast = Fmin;
        i10last = i10;
        F10last = F(i10);
    else
        relerr = -Inf;
        i10 = i10last;
        F10 = F10last;
    end
    dispf('... end in %0.4g s:  convergence rate = %0.4g',etime(clock,t),relerr)
end
[F10,iF10] = sort(F10);
i10 = i10(iF10);
n10 = length(F10);
stat = repmat(struct('int',[],'r2',[],'F',[],'p',[],'alpha',[],'x0',[],'h',[],'weight',[],'res',[],'resint',[],'iok',[]),1,n10);
b = zeros(2,n10);
for i=1:n10
    weight = w_(X0(i10(i)),H(i10(i)));
    iok = weight>0.01;
    [b(:,i),bint,r,rint,stats] = regress(y(iok),[o(iok) x(iok)],alpha);
%     [b(:,i),bint,r,rint,stats] = regress(weight.*y,[weight.*o weight.*x],alpha);
    stat(i) = struct(...
        'int',bint,...
        'r2',stats(1),...
        'F',stats(2),...
        'p',stats(3),...
        'alpha',alpha,...
        'x0',X0(i10(i))*sx+xm,...
        'h',H(i10(i))*sx,...
        'weight',weight,...
        'res',r,...
        'resint',rint,...
        'iok',i0ok(iok) ...
        );
end

% outputs
if ~nargout || ploton
    iplot = find([stat.r2]>0.6);
    nplot = min(maxplot,length(iplot));
    figure
    plot(x,y,'k-','linewidth',2)
    hold on
    col = jet(nplot);
    xplot = linspace(x(1),x(end),10)';
    oplot = ones(size(xplot));
    for i=1:nplot
        plot(x,[o x]*b(:,iplot(i)),'-','color',col(i,:))
        text(xplot,[oplot xplot]*b(:,iplot(i)),num2str(iplot(i)),'horizontalalignment','center','verticalalignment','middle','fontsize',10)
    end
end
if nargout, bout = b; end
if nargout>1, statout = stat; end

% Nested functions
    function Fval = F_(x0,h)
        % F_ = @(x0,h) (var(yn-[o xn]*b_(x0,h),w_(x0,nspline*h))/vy) ;
        weight  = w_(x0,nspline*h);
        nweight = length(find(weight>0.1));
        vy      = var(yn,weight);
        Fval =  (nweight+1)*(var(yn-[o xn]*b_(x0,h),weight)/vy);
        function coeff = b_(x0,h)
            % b_ = @(x0,h) regress(w_(x0,h).*yn,[w_(x0,h).*o w_(x0,h).*xn]);
            weight = w_(x0,h);
            coeff = regress(weight.*yn,[weight.*o weight.*xn]);
        end
        function weight = w_(x0,h)
            % w_ = @(x0,h) max(1-abs((xn-x0)/h).^3,0).^3;
            weight = max(1-abs((xn-x0)/h).^3,0).^3;
        end
    end
    function weight = w_(x0,h)
        % w_ = @(x0,h) max(1-abs((xn-x0)/h).^3,0).^3;
        weight = max(1-abs((xn-x0)/h).^3,0).^3;
    end

end
