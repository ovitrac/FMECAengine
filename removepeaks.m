 function [ynp,models] = removepeaks(x,y,segments,f,smooth,xtrials)
%  REMOVEPEAKS based on a regularized 3rd/4th degree interpolant polynomial with continuous values and 1st derivatives
%  Syntax: ynp = removepeaks(x,y,segments [,f,smooth])
%      x,y: data points read columnwise
% segments: [start1 stop1; start2 stop2...]
%        f: four degree coefficient (default = 0)
%   smooth: smooth parameter based on filtzero
%      ynp: output with peak removed
%  xtrials: try values x-xtrials
%
%
% Mathematical demonstration
%          syms a b c d x x1 x2 y1 yy1 y2 yy2
%         y = a *x^3+b*x^2+c*x+d;
%         y1exp = subs(y,x,x1);
%         yy1exp = subs(diff(y,x),x1);
%         y2exp = subs(y,x,x2);
%         yy2exp = subs(diff(y,x),x2);
%         res = solve(y1exp-y1,yy1exp-yy1,y2exp-y2,yy2exp-yy2,a,b,c,d);
%         p = simplify(subs(y,{a b c d},{res.a res.b res.c res.d}));
%         pretty(p)
%         pmatlab = matlabFunction(p);
%         %% degree 4 instead
%         syms f
%         y = f *x^4 + a *x^3+b*x^2+c*x+d;
%         y1exp = subs(y,x,x1);
%         yy1exp = subs(diff(y,x),x1);
%         y2exp = subs(y,x,x2);
%         yy2exp = subs(diff(y,x),x2);
%         res = solve(y1exp-y1,yy1exp-yy1,y2exp-y2,yy2exp-yy2,a,b,c,d);
%         p = simplify(subs(y,{a b c d},{res.a res.b res.c res.d}));
%         pretty(p)
%         pmatlab = matlabFunction(p);
%
%   See also: filtzero, ndf

%   MS 2.1 - 15/05/2012 - INRA\Olivier Vitrac - Xiaoyi Fang - rev. 27/06/12

% Revision history
% 27/06/12 add xtrials

% Definitions
fdefault = 0;
smooth_default = 5;
interpolant =  @(f,x,x1,x2,y1,y2,yy1,yy2)y2+x.*yy2-x2.*yy2+f.*x.^4-f.*x.^3.*x1.*2.0-f.*x.^3.*x2.*2.0-(x-x2).^3.*1.0./(x1-x2).^3.*(y1-y2).*2.0-((x-x2).^2.*(yy1+yy2.*2.0))./(x1-x2)+f.*x.^2.*x1.^2+f.*x.^2.*x2.^2+f.*x1.^2.*x2.^2+(x-x2).^2.*1.0./(x1-x2).^2.*(y1.*3.0-y2.*3.0+x.*yy1+x.*yy2-x2.*yy1-x2.*yy2)-f.*x.*x1.*x2.^2.*2.0-f.*x.*x1.^2.*x2.*2.0+f.*x.^2.*x1.*x2.*4.0;
xtrials_default = 0;

% Arg check
if nargin<3, error('3 inputs are required, ynp= removepeaks(x,y,segments,f,smooth)'), end
if nargin<4, f = []; end
if nargin<5, smooth = []; end
if nargin<6, xtrials = []; end
if isempty(f), f = fdefault; end
if isempty(smooth), smooth = smooth_default; end
if isempty(xtrials), xtrials = xtrials_default; end
[my,ny] = size(y);
if numel(x)>length(x), error('x must be a vector'), end
if length(x)~=my, error('x and y must have the same number of rows'), end
[nsegments,nends] = size(segments);
if nsegments==2 && nends==1, segments = segments'; nsegments = 1; nends = 2; end
if nends~=2, error('segments must be a nx2 array'), end
ntrials = length(xtrials);

% filtering and derivatives
yf = filtzero(y,smooth);   % filtered signal
dyfdx = ndf(x,yf,1);

% Correction for all segments
ynp = yf;
good = true(my,ntrials);
bad = false(my,ntrials);
if nargout>1, models = cell(nsegments,ny); end
for i=1:nsegments
    % find reference values
    bad0 = (x>=(segments(i,1))) & (x<=segments(i,2));
    for k=1:ntrials
        bad(:,k) = (x>=(segments(i,1)+xtrials(k))) & (x<=segments(i,2));
        good(:,k) = good(:,k) & ~bad(:,k);
    end
    % column-wise
    for j=1:ny
        tmp = cell(ntrials,1);
        crit = zeros(ntrials,1);
        %base = polyval(polyfit(segments(i,:),interp1(x(good(:,k)),yf(good(:,k),j),segments(i,:)),1),x(bad0));
        for k = 1:ntrials
            tmp{k} =  interpolant(f,...
                x(bad0),... x
                segments(i,1),segments(i,2),...
                interp1(x(good(:,k)),yf(good(:,k),j),segments(i,1)),... y1
                interp1(x(good(:,k)),yf(good(:,k),j),segments(i,2)),... y2
                interp1(x(good(:,k)),dyfdx(good(:,k),j),segments(i,1)),... yy1
                interp1(x(good(:,k)),dyfdx(good(:,k),j),segments(i,2)) ... yy2
                );
           crit(k) = sum(ndf(x(bad0),tmp{k},2).^2); %min(tmp{k}-base); % distance to be maximized
        end
        [~,k] = min(crit); %max(crit);
        if nargout>1
            xa = segments(i,1);
            xb = segments(i,2);
            y1 = interp1(x(good(:,k)),yf(good(:,k),j),segments(i,1));
            y2 = interp1(x(good(:,k)),yf(good(:,k),j),segments(i,2));
            yy1 = interp1(x(good(:,k)),dyfdx(good(:,k),j),segments(i,1));
            yy2 = interp1(x(good(:,k)),dyfdx(good(:,k),j),segments(i,2));
            models{i,j} = @(x)  interpolant(f,x,xa,xb,y1,y2,yy1,yy2);
        end
        ynp(bad0,j) = tmp{k};
        if xtrials(k)~=0, dispf('REMOVEPEAKS +%d',k),end
    end
end