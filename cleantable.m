function Tclean=cleantable(T,emptystring)
%CLEANTABLE cleans table/structure to be used BYKEYWORDS
%   Syntax: Tclean = cleantable(T [,emptystring))
%           Tclean=cleantable(T,emptystring)
%   Inputs:
%           T = structure created by LOADODS or XLSTBLREAD or table
%           emptystring: string to replace NaN values in string cell vectors (default='')
%
% See also: XLSTBLREAD, LOADODS, LOADODSPREFETCH, STRUCT2STRUCTTAB, STRUCTTAB2STRUCT, SUBSTRUCTARRAY, BYKEYWORDS

% MS 2.1 - 15/11/2015 - INRA\Olivier Vitrac - rev. 

% arg check
if nargin<1, error('one argument is required'), end
if nargin<2, emptystring = ''; end
if isstruct(T)
   if length(T)>1
       error('T must be a scalar structure, use STRUCTTAB2STRUCT if required to convert your strucutre')
   end
   T0 = T;
   isatable = false;
elseif istable(T)
    T0tab = table2struct(T);
    T0 = structtab2struct(T0tab);
    isatable = true;
else
    error('T must be a structure or a table')
end
n = unique(structfun(@length,T0));
if length(n)~=1, error('T fields must be array or numeric vectors'), end

%% fix bad columns and identify goodrows
badcolumns = structfun(@iscell,T0) & ~structfun(@iscellstr,T0);
fields = fieldnames(T0)';
goodrows = true(n,1);
for f = fields(badcolumns)
    tobeconverted = cellfun(@(x) isnumeric(x) && isnan(x) ,T0.(f{1}));
    if isatable
        T.(f{1})(tobeconverted) = {emptystring};
    else
        T0.(f{1})(tobeconverted) = {emptystring};
    end
    goodrows = goodrows & (cellfun(@ischar,T0.(f{1})) | tobeconverted);
end

%% extract good rows
if isatable
    Tclean = T(goodrows,:);
else
    T0tab = struct2structtab(T0);
    Tclean = structtab2struct(T0tab(goodrows));
end