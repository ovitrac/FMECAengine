function setprop(h,varargin)
%SETPROP retrieve properties from any object
%   Syntax: setprop(h,'property1',value1,'property2',value2...)
%           setfield(h,v)
%       h: nx1 or 1xn vector of handles
%       'property1', 'property2'...: valid properties for h
%       'value1', 'value2'...: corresponding values
%       v: structure with fields
%           v.property1 = value 1
%           v.property2 = value 2
%         >> array of structure generate an error


% MS 2.0 - 19/08/07 - Olivier Vitrac - rev.

% Revision history

% arg check
if nargin<2, error('syntax: setprop(h,''property1'',''value1'',''property2'',''value2'',...) or setfield(h,v)'), end
if ~all(ishandle(h)), error('some handles are invalid'), end
n = numel(h);
if isstruct(varargin{1})
    if nargin>2, error('with v0 structure, the syntax is v=getfield(h,v0)'), end
    prop = varargin{1};
    if length(prop)>1, error('v must be a structure and not an array of structure'), end
    proplist = fieldnames(prop)';
    m = length(proplist);
    varargin = cell(1,2*m);
    varargin(1:2:2*m-1) = proplist;
    varargin(2:2:2*m) = struct2cell(prop);
end
m = length(varargin);
if mod(m,2), error('arguments must be provided by pair properties/values'), end

% set
for i=1:n
    try
        set(h(i),varargin{:})
    catch
        for j=1:2:m-1
            try, set(h(i),varargin{j},varargin{j+1}), end
        end
    end
end
