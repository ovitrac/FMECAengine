function Tr = rentable(T,namexpr,repexpr)
%RENTABLE rename table/structure names according to a regular expression
% Example:
%{
    T = struct2table(struct('A',{'a' 'aa' 'aaa'},'B','b','C','c','Adup01','abis','Bbisdup01','bbis','Cbisdup02','cbis'));
    Tr = rentable(T,{'dup01$' 'dup02$'},{'' ''})
%}

% MS 2.1 - 03/08/2017 - INRA\Olivier Vitrac, FDA\Danielle Larese - rev. 

%arg check
if nargin<3, error('3 arguments are required'), end
isatable = istable(T);
if ~isstruct(T) && ~isatable, error('T must be a structure or a table'), end
if isatable
    isfield2 = @(T,f) ismember(f,T.Properties.VariableNames);
else
    isfield2 = @(S,f) isfield(S,f);
end

% for each field
Tr = T;
for f = fieldnames(T)'
    fr = regexprep(f{1},namexpr,repexpr);
    if ~isfield2(Tr,fr)
        Tr = renfield(Tr,f,fr);
    end
end