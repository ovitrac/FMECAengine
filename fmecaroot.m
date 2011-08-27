function root = fmecaroot(fmecadb,prop)
%FMECAROOT returns the root (according to prop) of each node (to be created with FMECAENGINE)
%  syntax: root = fmecaroot(fmecadb [,prop])
%       fmecadb: output object of FMECAengine
%          prop: 'parent' (default) or 'inherit'
%
%   See also: FMECAENGINE, FMECASINGLE, BUILDMARKOV, LOADFMECAENGINEDB, KEY2KEY

%MIGRATION 2.1 - 26/08/11 - INRA\Olivier Vitrac 

% arg check
if nargin<2, prop=''; end
if isempty(prop), prop = 'parent'; end

%scan
nodes = fieldnames(fmecadb)'; nnodes = length(nodes);
root = struct([]); % root = cell(nnodes,1);
for i=1:nnodes
    p = nodes{i};
    while ~isempty(fmecadb.(p).(prop))
        p = fmecadb.(p).(prop);
    end
    root(1).(nodes{i}) = p; %root{i} = p;
end