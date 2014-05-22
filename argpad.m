function y = argpad(x,n,mode,padvalue)
%ARGPAD pad elements in a vector/cell
%   y = argpad(x,n [,mode,padvalue])
%   x: numeric or cell array (char types are converted into {x})
%   n: required length
%   mode = keyword or padvalue
%       'last' (use the last value, default)
%       'periodic'/'first' (use the first value)

% MS 2.1 - 22/05/14 - INRA\Olivier Vitrac - rev.


% default
modedefault = 'last';

% arg check
if nargin<2, error('2 arguments are required'), end
if nargin<3, mode = ''; end
if isempty(x), y = x; return, end
if ischar(x), x = {x}; end
if isempty(mode), mode = modedefault; end
if ischar(mode)
    if strcmpi(mode,'last'), padvalue = x(end);
    elseif strcmpi(mode,'first') || strcmpi(mode,'periodic'), padvalue = x(1);
    else error('''%s'' is an unknown mode (available modes: ''last'' or ''first''/''periodic'')',mode)
    end
end
if iscell(x) && ~iscell(padvalue), padvalue = {padvalue}; end
m = numel(x);
if length(x)<m, x = x(:); end
if n<=m, y=x; return, end
if size(x,2)==1
    y = [x;padvalue(ones(n-m,1))];
else
    y = [x,padvalue(ones(1,n-m))];
end
