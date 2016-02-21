function combout = createchoosek(entity,nchoose,nosymmetry)
%CREATECHOOSEK creates a matrix of all combinations of entities
%   Syntax: comb = createchoosek(entity,nchoose [,nosymmetry])
%         entity = raw array of entities (numbers or letters)
%        nchoose = number of entities to pick
%     nossymetry = flag to remove symmetric permutations
%
%   Examples:
%     createchoosek('CE',5)
%     createchoosek('CE',5,true)
%     createchoosek([0 1],5,true)
%     createchoosek([0 1e6],5,true)
%     createchoosek([0 NaN],5,true)
%
% Example for the paper PVAc (15/01/2015, repetition of Table 5.13 of page 287 in Mai's Thesis)
%{
    monomers = 'AE';
    reference = 1;
    for n=1:6
        dispf('n=%d',n)
        s = createchoosek(monomers,n,true);
        y = cellfun(@(p) length(find(p==monomers(reference))),s);
        for yu = unique(y)'
            scurrent = s(y==yu); nu = length(scurrent);
            [num,den] = rat(yu/n); if den==1, ytxt = sprintf('%d\t',num); else ytxt = sprintf('%d/%d',num,den); end
            %txt = cellfun(@(i,t) sprintf('%s:%s, ',i,t), arrayfun(@(i) char(96+i),1:nu,'UniformOutput',false)',scurrent,'UniformOutput',false);
            txt = cellfun(@(i,t) sprintf('%s, ',t), arrayfun(@(i) char(96+i),1:nu,'UniformOutput',false)',scurrent,'UniformOutput',false);
            dispf('\ty=%s\t%s (%d)',ytxt,regexprep([txt{:}],'\,\s$',''),nu)
        end
    end
%}

%MS 2.1 - 07/07/14 - INRA\Olivier Vitrac, Wafa Guiga, Mai Nguyen - rev. 15/01/2015


% revision History
% 15/01/2015 add advanced example

% default
nchoose_default = 1;
nosymmetry_default = false;
alphabet = char(32:127); % keep them continuous
nalphabet=length(alphabet);

% argcheck
if nargin < 1, error('At least one argument is required'), end
if nargin < 2, nchoose = []; end
if nargin < 3, nosymmetry = []; end
if isempty(nchoose), nchoose = nchoose_default; end
if isempty(nosymmetry), nosymmetry = nosymmetry_default; end
entity = entity(:)'; nentity = length(entity);
ischarentity = ischar(entity);

% recursion if nosymmetry is combined with numeric entity
if ~ischarentity && nosymmetry
   if nentity>nalphabet, error('the maximum numb of entities (%d) is exceeded',nentity), end
   [entityu,ientityu] = unique(entity);
    combout = entityu(char(createchoosek(alphabet(ientityu),nchoose,nosymmetry))-alphabet(1)+1);
    return
end

% do the job
entitieswithrepetitions = kron(ones(1,nchoose),entity);
if ischarentity, entitieswithrepetitions = char(entitieswithrepetitions); end
combout = entitieswithrepetitions(nchoosek(1:nentity*nchoose,nchoose));
if ischarentity
    combout=unique(cellstr(combout));
else
    combout=unique(combout,'rows');
end

% remove symmetry if required
if nosymmetry && ischarentity
    combout = unique(removesymmetryinstrings(combout));
end