function X0=renfield(X,f,f0)
%RENFIELD rename a field in a structure (can be applied to several fields)
%   syntax1: X0=renfield(X,f,f0)
%       X: a structure or an array of structures
%       f, f0: fields names (string or cell arrays of string)
%       X0: structure with fields in the similar order as X
%   syntax2: X0=renfield(X,f0)
%       It assumes f = fieldnames(X) and length(f0)=length(f)

%   example 1:
%   X = struct('a',{'a1' 'a2'},'b',{'b1' {'b21' 'b22'}},'c','...')
%   X0 = renfield(X,{'a' 'c'},{'A' 'C'})
%
%   example 2: as above but convert all fieldnames into UPPERCASE
%   X0 = renfield(X,upper(fieldnames(X)))
%
%   SEE ALSO: REPLACESYNONYMS, BYKEYWORDS, XLSTBLREAD, LAMDUMPREAD (any functon involving STRUCT)
%
%   Dependencies for distribution: NONE

% MS 2.1 - 23/01/08 - INRA\Olivier Vitrac rev. 03/08/17

% revision history
% 30/01/08 remplace length(X) by numel(X)
% 03/11/16 fix nargchk (for future Matlab versions)
% 03/08/2017 accept tables

% arg check
if nargin<2 || nargin>3, error('syntax: X0=renfield(X,f,f0) or X0=renfield(X,f0)'), end
isatable = istable(X);
if ~isstruct(X) && ~isatable, error('X must be a structure or a table'), end
if nargin==2, f0=f; f=fieldnames(X); end
if ~iscell(f), f = {f}; end
if ~iscell(f0), f0 = {f0}; end
if ~iscellstr(f), error('f must be a string or cell array of strings'), end
if ~iscellstr(f0), error('f0 must be a string or cell array of strings'), end
if numel(f)~=numel(f0), error('f and f0 must be of a same size'), end

% scan
if isatable
    fX = X.Properties.VariableNames;
    m = 1;
else
    fX = fieldnames(X);
    m = numel(X); % for array of structure
end
for i=1:length(fX)
    j = find(ismember(f,fX{i}));
    if any(j)
        if length(j)>1, error('the field ''%s'' is duplicated in f',fX{i}), end
        [X0(1:m).(f0{j})] = deal(X.(fX{i}));
    else
        [X0(1:m).(fX{i})] = deal(X.(fX{i}));
    end
end

if isatable, X0 = struct2table(X0); end