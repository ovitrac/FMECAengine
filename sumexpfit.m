function R=sumexpfit(t,y,ordre,Cinf,options)
% SUMEXPFIT calcule le développement à l'ordre n de y(t) en somme d'exponentielles décroissantes/croissantes
%		R=sumexpfit(t,y,ordre,Cinf,options)	
%
%       See also SUMEXPDER, SUMEXPDER


% Woodox 1.21 - 02/04/01 - Olivier Vitrac - rev. 01/07/13

% revision history
% 28/04/06 fit y(end)-y if y is increasing (i.e. switch type from 'decexp' to 'incexp')
% 02/10/09 fix lambda(1) guess with -regress(log(y),t)
% 05/11/12 fix sorting with negative amplitudes
% 01/07/13 drop Diagnostics, not anymore supported

% Inputs control
warnoff = 1;
if nargin<3, ordre = 2; end
if nargin<4, Cinf = 0; end
if nargin<5, options = optimset; end
% else warnoff = strcmp(optimget(options,'Diagnostics'),'off'); end
if isempty(Cinf), Cinf = 0; end
t = t(:); y = y(:);
if y(end)>y(1)
    type = 'incexp';
    Cend = y(end);
    y = y(end)-y;
else
    type = 'decexp';
    Cend = 0;
end
t0 = min(t); tend = max(t); yend = mean(y(t==tend));
if yend>.1*(max(y)-min(y)); Cinf = yend; end
if warnoff, warning off, end %#ok<WNOFF>

% First guess
lambda=1./((1:ordre)*(tend-t0));
simple = [ones(size(t)) t];
if all(y>0)
    logy =log(y);
    guess = regress(logy,simple);
    lambda(1) = max(-guess(2),-regress(logy,t));
    for i=2:ordre
        logy =logy-(guess(1)+guess(2)*t);
        guess = regress(logy,simple);
        lambda(i) = max(-regress(logy,t),max(lambda(i),-guess(2)));
    end
end
%[err,C]=sumexperr(lambda,t,y,Cinf);
% Solving
[lambda,~,flag] = fminsearch(@sumexperr,lambda,options,t,y,Cinf);
% Solution
[err,C]=sumexperr(lambda,t,y,Cinf);
lambda(C==0) = NaN;

% Outputs control
[~,ind] = sort(abs(C),'descend');
R = struct( 'type', type,...
            'order', length(find(abs(C)>sqrt(eps))),...
            'C',     C(ind)',...
            'l',    lambda(ind),...
            'Cinf', Cinf,...
            'Cend', Cend,... used if type = incexp as ytrue = Cend-ypred
            'err',  err,...
            'r2',   1-err/(var(y)*(length(y)-1)),...
            'flag', flag );
if warnoff, warning on, end %#ok<WNON>



function [err,C] = sumexperr(lambda,t,y,Cinf)
% non linear part
m = length(lambda);
LAMBDAS = zeros(length(t),length(lambda));
for i = 1:m, LAMBDAS(:,i) = exp(-lambda(i)*t); end
% rank check
M = LAMBDAS'*LAMBDAS;
if ~any(isinf(M(:))) && ~any(isnan(M(:)))
    r = rank(M,max(size(M))*norm(M)*sqrt(eps));
    if r<m, LAMBDAS = LAMBDAS(:,1:r); end
else
    r = Inf;
end
C = LAMBDAS\(y-Cinf);
ye = LAMBDAS*C;
err = norm(ye-(y-Cinf));
if r<m, C(r+1:m) = 0; end