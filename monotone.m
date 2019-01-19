function [pos,width,height] = monotone(X,type,zero,keyword)
% MONOTONE looks for monotonous segments in vector X
%		ex. iX = monotone(X)
%		options : [pos,width,height] = monotonie(X [,type,zero,keyword])
%			inputs :
%				X = column vector
%				type = '+' (default), '-' or '0' ('+-' is also implemented)
%				zero = maximum value of dX (diff(X)) associated to 0)
%            keyword =  '' (default) - as above
%                      'struct' to retrieve all output arguments in a structure
%                      'table' returns a table instead of a struct
%			outputs :
%				pos = beginning position of the monotonous segment (index)
%				width = width of the segment (index unit)
%				amplitude/height of the segment
%               when struct is used: struct('position',pos,'width',width,'height',height)

% Thermique 1.0 - 30/04/01 - Olivier Vitrac (source : WoodOx 1.22) - rev. 07/08/18

% Revision History
% 26/07/10 add '+-'
% 28/09/12 change default zero from 1e5*eps to max(10*eps,(max(X)-min(X))/1e6)
% 23/05/13 accept empty zero
% 07/08/18 sort positions if '+-' is used, English translation, add keyword
% 09/08/18 fix 'struct' and 'table' when '+' or '-' is used alone

if nargin<2, type='+'; end
if nargin<3, zero = []; end
if nargin<4, keyword = ''; end
if isempty(zero), zero = max(10*eps,(max(X)-min(X))/1e6); end
if strcmp(type,'+-') || strcmp(type,'-+')
    [pp,lp,ap] = monotone(X,'+',zero);
    [pm,lm,am] = monotone(X,'-',zero);
    %pos = [pp;pm]; % old implementation
    % new implementation (note sort and isort)
    [pos,isort] = sort([pp;pm],'ascend');
    if ischar(keyword) && (strcmpi(keyword,'struct') || strcmpi(keyword,'table'))
        largtmp = [lp;lm];
        amptmp = [ap;am];
        pos = struct('position',pos,'width',largtmp(isort),'height',amptmp(isort));
        if strcmpi(keyword,'table'), pos = struct2table(pos); end
    else
        if nargout>1, width = [lp;lm]; width = width(isort); end
        if nargout>2, height = [ap;am]; height = height(isort); end
    end
    return
end

X = X(:);
dX=diff(X);
S = [(dX>zero)-(dX<-zero)];
switch upper(type)
case '+', indS = find(S>0);
case '-', indS = find(S<0);
case '0', indS = find(S==0);
end
%dindS = [2;diff(indS)]; % not used 28/09/12
if any(indS)
	il = indS([2;diff(indS)]>1);
	ir = indS([diff(indS);2]>1)+1;
else
   il = []; ir = [];
end

if ~nargout && any(il)
   hold on
   plot(X,'b-'), plot(il,X(il),'ro',ir,X(ir),'ms')
   stem(il,X(ir)-X(il),'gd')
else
   if  nargout>0, pos = il; end
   if  nargout>1, width = ir-il +1; end
   if  nargout>2, height = X(ir)-X(il); end
   if ischar(keyword) && (strcmpi(keyword,'struct') || strcmpi(keyword,'table'))
       pos = struct('position',il,'width',ir-il+1,'height',X(ir)-X(il));
       if strcmpi(keyword,'table'), pos = struct2table(pos); end
   end
end