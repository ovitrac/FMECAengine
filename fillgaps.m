function filled = fillgaps(t,y,varargin)
%FILLGAPS fill gaps in kinetics based on a minimum sampling time (a power law correction is applied, default exponent matches diffusive behavior)
%   filled = fillgaps(t,y,'property',value)
%
%   List of property/value (default)
%             dt: minimum time interval between the first and the second point (max(0,min(diff(t)))/2)
%       exponent: exponent value to the variable transform (t.^exponent), (default=0.5)
%            tol: tolerance for spaps (default = 1e-4; i.e. 0.01% of the variance)
%
%   Output: filled.t, filled.y

% Migration 2.0 - 04/01/2019 - INRA\Olivier Vitrac - rev.

% Default
default = struct(...
    'dt',[], ... same unit as t
    'exponent',0.5,...
    'tol',1e-4 ...
    );

% argcheck
if nargin<2, error('2 arguments are required: yfilled = fillgaps(t,y,''property'',value)'), end
[t,ind] = sort(t(:)); y = y(ind); n = length(t);
if length(y)~=n, error('t and y should have the same length'), end
o = argcheck(varargin,default);
if isempty(o.dt), o.dt = max(0,min(diff(t)))/2; end

% do the job
ttransform = t.^o.exponent;
dttransform = o.dt^o.exponent;
igap = find(diff(ttransform)>dttransform);
if ~isempty(igap)
    t2add = arrayfun(@(t1,t2) (t1+dttransform:dttransform:t2-dttransform)',ttransform(igap),ttransform(igap+1),'UniformOutput',false);
    t2add = cat(1,t2add{:});
    sp = spaps(ttransform,y,o.tol*var(y)*length(y));
    y2add = fnval(sp,t2add);
    yall = [y;y2add];
    [tall,ind] = unique([ttransform;t2add]);
    filled = struct('t',tall.^(1/(o.exponent)),'y',yall(ind));
else
    filled = struct('t',t,'y',y);
end