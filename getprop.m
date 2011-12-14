function v=getprop(h,varargin)
%GETPROP retrieve properties from any object
%   Syntax: v=getprop(h,'property1','property2'...)
%           v=getprop(h,{'property1','property2'...})
%           v=getfield(h,v0)
%           v=getprop(h)
%       h: nx1 or 1xn vector of handles
%       'property1', 'property2'...: valid properties for h
%       v, v0: structure with fields
%           v.property1 = value 1
%           v.property2 = value 2
%       >> readonly properties and handles properties are excluded by default (cannot be copied with setprop)
%          use v=getprop(h,'-override',[...]) to override this behavior
%       >> non existing property generate an empty output

% MS 2.0 - 19/08/07 - Olivier Vitrac - rev. 04/09/07

% Revision history
% 03/09/07 add v=getprop(h), remove readonly properties
% 04/09/07 add '-override', add getprop(handle,{'prop1' ... 'propn'}), try

% readonly, incompatible, handles properties (all properties wich cannot be copied with setprop)
readonly = {'beingdeleted' 'children' 'parent' 'tightinset' 'title' 'type' 'xlabel' 'ylabel' 'zlabel' ...
'currentpoint'  'currentaxes' 'currentcharacter' 'currentmenu' 'currentobject' 'uicontextmenu' ...
'FixedColors' 'location' 'extent' ...
};
override = '-override';

% arg check
if nargin<1, error('syntax: v=getprop(h,''property1'',''property2''...) or v=getfield(h,v0)'), end
if ~all(ishandle(h)), error('some handles are invalid'), end
effnargin = nargin;
if nargin>1 && ischar(varargin{1}) && strcmp(varargin{1},override)
    effnargin = effnargin-1;
    varargin = varargin(2:end);
    readonly = {};
end
if effnargin<2,
    s=get(h); f=fieldnames(s); va = struct2cell(s);
    [f,ind]=setdiff(lower(f),readonly);
    v = cell2struct(va(ind,:),f,1);
    return
end
if effnargin==2 && iscell(varargin{1}), varargin = varargin{1}; end
n = numel(h);
if isstruct(varargin{1})
    if nargin>2, error('with v0 structure, the syntax is v=getfield(h,v0)'), end
    varargin = fieldnames(varargin{1})';
end
m = length(varargin);

% get
v = repmat(cell2struct(repmat({[]},1,m),varargin,2),n,1);
for f=setdiff(lower(varargin),readonly)
    for i=1:n
        try, v(i).(f{1})=get(h(i),f{1}); end
    end
end