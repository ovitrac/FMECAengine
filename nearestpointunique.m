function ind=nearestpointunique(val,list,mode,tolerant)
%NEARESTPOINTUNIQUE find unique indices in a list matching "in some way" a set of values
%   syntax: ind=nearespointunique(val,list [,mode,tolerant])
%       val = values to search
%      list = reference list
%      mode = 'nearest', 'near', 'absolute' (default = 'nearest')
%  tolerant = force nearestpoint instead of nearestpointunique if val is larger than list (rendundancy is introduced) 
%       ind = index of matching values in list 
%
%   See also: nearestpoint
%
%   Example: nearestpointunique([4 5.1 4.1 4.2 4.4 5],1:10)
%   gives [4     6     3     2     7     5]

% MS 2.1 - 06/12/2015 - INRA\Olivier Vitrac - rev. 02/12/2017

% default
mode_default = 'nearest';

% arg check
if nargin<2, error('2 arguments are required'); end
if nargin<3, mode = ''; end
if nargin<4, tolerant = false; end
if isempty(mode), mode = mode_default; end
m = numel(val); n=numel(list);
if m>n
    if ~tolerant
        error('the number of elements in val (%d) must be smaller or equal than the number of elements in list (%d)',m,n)
    else
        warning('the number of elements in val (%d) must be smaller or equal than the number of elements in list (%d)',m,n)
        ind=nearestpoint(val,list,mode);
        return
    end
end

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