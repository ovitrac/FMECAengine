function f = fmecasortdb(fmecadb)
%FMECASORTDB sorts a FMECAengine database in a similar way it was calculated
%   syntax: fmecadbsorted=fmecasortdb(fmecadb)

% Migration 2.1 (Fmecaengine v0.497) - 10/04/2011 - INRA\Olivier Vitrac - rev. 

% arg check
if ~isstruct(fmecadb), error('the argument must be a structure'), end
steps =  fieldnames(fmecadb); nsteps = length(steps);
if ~all(cellfun(@(x) isfield(fmecadb.(x),'Bi'),steps)), error('the database must be created with FMECAengine'), end

% paths
[~,paths] = buildmarkov(fmecadb);
[~,pathorder] = sort(cellfun(@(x) length(x),paths),'descend');
paths = paths(pathorder);
orderedlist = cat(2,paths{:});

% reorder
% table(i) gives the order of step i
% table2(order) = i (reverse table)
[table,table2] = deal(zeros(nsteps,1));
rank = 0;
for i=orderedlist
    if table(i)==0, rank=rank+1; table(i)=rank; end
end
table2(table)=1:nsteps; % bijection

% output
f = orderfields(fmecadb,table2);
