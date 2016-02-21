function ind=nearestpointunique(val,list,mode)
%NEARESTPOINTUNIQUE find unique indices in a list matching "in some way" a set of values
%   syntax: ind=nearespointunique(val,list [,mode])
%       val = values to search
%      list = reference list
%      mode = 'nearest', 'near', 'absolute' (default = 'nearest')
%       ind = index of matching values in list 
%
%   See also: nearestpoint
%
%   Example: nearestpointunique([4 5.1 4.1 4.2 4.4 5],1:10)
%   gives [4     6     3     2     7     5]

% MS 2.1 - 06/12/2015 - INRA\Olivier Vitrac - rev. 

% default
mode_default = 'nearest';

% arg check
if nargin<2, error('2 arguments are required'); end
if nargin<3, mode = ''; end
if isempty(mode), mode = mode_default; end
m = numel(val); n=numel(list);
if m>n, error('the number of elements in val (%d) must be smaller or equal than the number of elements in list (%d)',m,n), end

% guess the best order
crit = abs(val-list(nearestpoint(val,list,mode)));
[~,order] = sort(crit(:)','ascend'); %#ok<UDIM>

% do
ind = zeros(size(val));
ok = true(n,1);
for i=order
    jvalid = find(ok);  
    ind(i) = jvalid(nearestpoint(val(i),list(ok),mode));
    ok(ind(i))=false;
end