function s = structtab2struct(sarray,dim,NonUniformOutput)
%STRUCTTAB2STRUCT converts a struct array into a structure with array/cell field contents (reverse operation of struct2stucttab)
%   syntax: s = structtab2struct(sarray [,dim,NoNUniformOutput])
%           sarray: structure array
%              dim: dimension along which the collapse is applied
% NoNUniformOutput: flag, to be set to true when all fields do not have the same size
%
%   example: structtab2struct(dir)
%
%   example with dissimilar elements:
%{
    a.f1 = [1;2]; a.f2 = {{'a1';'a2';'a3'};'b'}; a.f3 = {(1:3)'; 2:3};
    b.f1 = [3;4;5]; b.f2 = {'e'; {'f';'g'}; {'h' 'i'}}; b.f3 = {3; (4:5)' ; (40:60)'};
    as = struct2structtab(a);
    bs = struct2structtab(b);
    ab = structtab2struct([as;bs],[],true);
%}
%
%   See also STRUCTTAB2STRUCT, SUBSTRUCTARRAY, CLEANTABLE, FLATENSTRUCTTAB, TAB2CSB


% MS 2.1 - 29/11/13 - INRA\Olivier Vitrac - rev. 24/11/2015

% Revision history
%  05/06/2014 add see also
%  09/03/2015 add NonUniformOutput
%  10/03/2015 remove uncell when NonUniformOutput=true
%  24/11/2015 update See also section

% Default
dim_default = 1;
NonUniformOutput_default = false;

% argcheck
if nargin<1, error('one argument is required'), end
if ~isstruct(sarray), error('the first argument must be a structure'), end
if nargin<2, dim = []; end
if nargin<3, NonUniformOutput = []; end
if isempty(dim), dim = dim_default; end
if isempty(NonUniformOutput), NonUniformOutput = NonUniformOutput_default; end


% non-uniform outputs
fn = fieldnames(sarray);
if NonUniformOutput % non uniform outputs (i.e. non-scalar outputs) need to protected
    nsarray = length(sarray);
    for j=1:length(fn)
        i=1; refsiz = size(sarray(i).(fn{j})); uniform = prod(refsiz)==1;
        while uniform && i<nsarray
            i=i+1; uniform = matcmp(size(sarray(i).(fn{j})),refsiz);
        end
        if ~uniform
            for i=1:nsarray, sarray(i).(fn{j}) = {sarray(i).(fn{j})}; end
        end
    end
    % collapsing
    s = cell2struct(cellfun(@(f) cat(dim,sarray.(f)),fn,'UniformOutput',false),fn);
else
    % collapsing
    s = cell2struct(cellfun(@(f) substructarray(sarray,f,dim),fn,'UniformOutput',false),fn);
end

