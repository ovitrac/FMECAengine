function t=struct2structtab(s)
%STRUCT2STRUCTTAB cast a structure into a structure array
%      t = struct2structtab(s)

% MS 2.1 - 20/01/12 - INRA\Olivier Vitrac rev.

% argcheck
if nargin<1 || nargin>1, error('one argumentis required'), end
if ~isstruct(s), error('the argument must be a structure'), end
if isempty(s) || numel(s)>1, error('structure arrays or empty structure are not authorized'), end
if any(diff(structfun(@numel,s))), error('all fields must be of a same size'), end
f = fieldnames(s);

% convert numeric fields
fn = f(structfun(@isnumeric,s));
for i=1:length(fn)
    s.(fn{i}) = num2cell(s.(fn{i})(:),2);
end
v = struct2cell(s);

% recast
tmp = [f v]';
t = struct(tmp{:});