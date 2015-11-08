function y = argpad(x,n,mode,padvalue,truncate)
%ARGPAD pad elements in a vector/cell
%   y = argpad(x,n [,mode,padvalue,truncateflag])
%   x: numeric or cell array (char types are converted into {x})
%   n: required length
%   mode = keyword or padvalue
%       'last' (use the last value, default)
%       'first' (use the first value)
%       'periodic' (periodic values)
%       'symmetric' (symmetric values)
%   truncateflag: forces truncation to n (default = true)

% MS 2.1 - 22/05/14 - INRA\Olivier Vitrac - rev. 07/11/2015

% Revision history
% 01/11/2015 add truncate (by default, truncation to n arguments is applied, it was not the case before)
% 01/11/2015 fix padvalue when supplied by the user
% 07/11:2015 full implementation of periodic and symmetric bounary conditions

% default
truncate_default = true;
modedefault = 'last';

% arg check
if nargin<2, error('2 arguments are required'), end
if n<0 || n~=round(n), error('n must be a positive integer'); end
if nargin<3, mode = ''; end
if nargin<4, padvalue = []; end
if nargin<5, truncate = []; end
if isempty(x), y = x; return, end
if ischar(x), x = {x}; end
if isempty(mode), mode = modedefault; end
if ~ischar(mode) || ~ismember(lower(mode),{'first','last','periodic','symmetric'}), error('mode must first/last/periodic/symmetric'), end
if isempty(truncate), truncate = truncate_default; end
m = numel(x); if length(x)<m, x = x(:); end
mpad = length(padvalue);
if mpad<1
    if (m>=n)
        padvalue = [];
    elseif strcmpi(mode,'last') % default
        padvalue = repmat(x(end),n-m,1);
    elseif strcmpi(mode,'first')
        padvalue = repmat(x(1),n-m,1);
    elseif strcmpi(mode,'periodic')
        padvalue = repmat(x(:),ceil((n-m)/m)*m,1);
    elseif strcmpi(mode,'symmetric')
        xtmp = x(:);
        padvalue = repmat([flipud(xtmp(1:end-1));xtmp(2:end)],ceil((n-m)/(2*(m-1)))*2*(m-1),1);
    end
else
    if (m>=n)
        padvalue = [];
    elseif strcmpi(mode,'last') % default
        padvalue = repmat(padvalue(end),n-m,1);
    elseif strcmpi(mode,'first')
        padvalue = repmat(padvalue(1),n-m,1);
    elseif strcmpi(mode,'periodic')
        padvalue = repmat(padvalue(:),ceil((n-m)/mpad)*mpad,1);
    elseif strcmpi(mode,'symmetric')
        if mpad<2, error('2 pad values are required at least to enable symmetry boundary conditions'); end
        xtmp = padvalue(:);
        padvalue = repmat([flipud(xtmp(1:end-1));xtmp(2:end)],ceil((n-m)/(2*(mpad-1)))*2*(mpad-1),1);
    end
    
end
padvalue = padvalue(1:n-m);

if iscell(x) && ~iscell(padvalue), padvalue = {padvalue}; end

% Do the job
if n<=m
    y=x;
else
    if size(x,2)==1
        y = [x;padvalue];
    else
        y = [x,padvalue'];
    end
end
if truncate, y = y(1:n); end