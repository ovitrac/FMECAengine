function [b,bstat] = xyrobustregress(x,y,alpha)
%ROBUSTREGRESS combination of regress(y(:),[ones(size(x,1),1) x(:)]) robustfit(x(:),y(:)) to generate confidence interval not sensitive to outliers
% Syntax: b = robustregress(x,y [,alpha])
%         [b,bstat] = robustregress(...)
%     x,y: nx1 vector
%   alpha: tp compute 100*(1-alpha)% confidence level (default = 0.05)
%       b: [intercept slope]
%   bstat: [ interceptmin interceptmax
%            slopemin     slopemax]

% MS 2.1 - 02/03/12 - Olivier Vitrac - rev. 

if nargin<2, error('2 arguments are required'), end
if nargin<3, alpha = 0.05; end
if numel(x) ~=numel(y), error('x and y must have the same size'); end
x = x(:); y = y(:); n = length(x);
[b1,b1stat] = regress(y,[ones(n,1) x]);
[b2,b2stat] = robustfit(x,y);
tval = tinv((1-alpha(1)/2),b2stat.dfe);
b2stat = [b2-tval*b2stat.se, b2+tval*b2stat.se];
performance = mean([(b1stat(:,2)-b2)./abs(b1) (b2stat(:,2)-b2)./abs(b2)]);
if performance(2)<performance(1)/1.05 % confidence interval reduced by at least 5% when robustfit is used instead regress
    b = b2;
    if nargout>1, bstat = b2stat; end
else
    b = b1;
    if nargout>1, bstat = b1stat; end
end
