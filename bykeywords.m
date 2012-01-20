function [X0,foundindexout]=bykeywords(X,f,varargin)
%BYKEYWORDS index structures or cells by keywords, the result is a similar but indexed object (useful with XLSTBLREAD)
%   X0 = bykeywords(X,key,keyvalues)
%   X0 = bykeywords(X,{key1 key2...},keyvaluesforkey1,keyvaluesforkey2,...)
%   [X0,foundindex] = bykeywords(...)
%          X: simple structure (i.e. X.fieldi = mx1 vector or mxn array) or mxn cell array
%             structure example: s=struct('kw',{{'kw1' 'kw2' 'kw3'}},'field1',1:3,'field2',10:10:30)
%             cell example1: c={{'kw1' 'kw2' 'kw3'} 1:3 10:10:30}
%             cell example2: c={'kw1' 'kw2' 'kw3'; 1 2 3; 10 20 30}'; % note TRANSPOSE to make a columwise table (default)
%              structure example: s=struct('kw1',{{'kw11' 'kw12' 'kw13'}},'kw2',{{'kw21' 'kw22' 'kw23'}},'field1',1:3,'field2',10:10:30)
%        key: field or column index (between 1 and n)
%  keyvalues: string or cell array of string
%         X0: index structure or cell array based on X and keys
% foundindex: found index
%
%     Duplicate keyvalue in X generates an error except if a single KEYVALUE is used
%     Redundancy in keyvalues generates an error except if a single INDEX is usedis accepted only for unique index.
%     Missing KEYVALUEs generate an error except if a single INDEX is used, it is accepted only for an unique index.
%     KEY and KEYVALUEs are case sensitive. If no match is found, a non-case sensitive comparison is tested
%
%   See also: XLSTBLREAD
%
%   Advanced examples from Guillaume thesis
%
%     Example1: extraction of multiple indices from a single combination of keyvalues: 'BHT'x'Ethanol'
%       tab=xlstblread('c:\data\olivier\e-mail\attach\Chi_database19','complete');
%       subtab = bykeywords(tab,{'Base' 'Screen'},'BHT','Ethanol');
%       >> result overview
%         Base: {2x1 cell}
%         Screen: {2x1 cell}
%       Energies: {2x1 cell}
%        Chi313K: [2x1 double]
%       Emix313K: [2x1 double] ... etc (Note that 2 entries are available for each field)
%
%     Example2: extraction of indices from a single key and several keyvalues
%       tab=xlstblread('c:\data\olivier\inra\codes\partition\dbK');
%       subtab = bykeywords(tab,'molecule',{'decane' 'undecane' 'dodecane' 'tridecane'})
%       >> result overview
%        molecule: {4x1 cell}
%        connolly: [4x1 double]
%     connollystd: [4x1 double]
%           molar: [4x1 double]
%         density: [4x1 double]
%              MW: [4x1 double]
%              nC: [4x1 double]
%
%   Dependencies for distribution: none
%
%   See also: XLSTBLREAD

% MS 2.0 - 20/01/08 - INRA\Olivier Vitrac - rev. 29/01/08

% Revision history
% 29/01/08 add finindex, new rules for missing KEYVALUES

% Definitions
reshapeon =false;

% arg check
if nargin<3, error('syntax: X0 = bykeywords(X,key,keyvalues)\nX0 = bykeywords(X,{key1 key2...},keyvaluesforkey1,keyvaluesforkey2,...)'), end
if isstruct(X)
    Xfields = fieldnames(X);
    if ~iscell(f), f = {f}; end
    ndim = length(f);
    f0 = zeros(1,ndim);
    for idim = 1:ndim
        if ~ischar(f{idim}), error('the key must be a string (single index) or cell array of strings (several indices)'), end
        j = find(ismember(Xfields,f{idim}));
        if isempty(j), error('the key ''%s'' is missing',f{idim}), end
        f0(idim) = j;
    end
    f = f0;
    X = struct2cell(X);
else
    Xfields = {};
    if  ~iscell(X), error('invalid object, only structures and cell arrays are accepted'), end
    if all(~cellfun(@iscell,X(:)))
       ncol = size(X,2);
       X0 = X; X = cell(1,ncol);
       for k=1:ncol, X{k} = X0(:,k); end % reshape
       reshapeon = true;
    end
    if rem(f,f) || f<=0, error('the index must be a positive integer'), end
end
ncol = numel(X);
if any(f>ncol), error('out of range index/indices'), end

% check keyvalues
keyvalues = varargin;
if length(keyvalues)~=ndim, error('the size of key indexing is inconsistent'); end
for idim = 1:ndim
    if ~iscell(keyvalues{idim}), keyvalues{idim} = {keyvalues{idim}}; end;
    if idim>1 && length(keyvalues{idim})~=length(keyvalues{1}), error('the length of keyvalues must be equal'), end
    if ~iscellstr(X{f(idim)}), error('the %dth dim must be a string or cell array of strings (mixed types are not permitted)',idim), end
    if any(cellfun('isempty',X{f(idim)})), error('some keyvalues along the %dth dim are empty.',idim), end
end
nkeyvalues = length(keyvalues{1});

% lookfor keyvalues, merge mutiple indices into a single one (inconsistencies generate errors)
i = zeros(nkeyvalues,1);
for j=1:nkeyvalues % for all vector of keyvalues
    ik = find(ismember(X{f(1)},keyvalues{1}(j)));
    if ndim>1
        for idim = 1:ndim
            tmp = find(ismember(X{f(idim)},keyvalues{idim}(j)));
            if ~any(tmp)
                tmp = find(ismember(lower(X{f(idim)}),lower(keyvalues{idim}(j))));
                if any(tmp), warning('a match for the KEYVALUE ''%s'' along DIM %d is found using a case insensitive search',keyvalues{idim}{j},idim), end
            end
            if ~any(tmp) && nkeyvalues>1, error('the KEYVALUE ''%s'' is missing along dim %d',keyvalues{idim}{j},idim), end
            ik = intersect(ik,tmp);
        end
    elseif ~any(ik)
        ik = find(ismember(lower(X{f(1)}),lower(keyvalues{1}(j)))); % case insensitive search
        if any(ik), warning('a match for the %dth KEYVALUE (''%s'') is found using a case insensitive search',j,keyvalues{1}{j}), end
    end
    if ~any(ik)
        warning('no match for the %dth KEYVALUE::',j)
        for idim = 1:ndim, disp(sprintf('dim %d: %s',idim,keyvalues{idim}{j})), end
    else
        if nkeyvalues>1
            if length(ik)>1, error('%dth index is inconsistent',j), end
            i(j) = ik;
        else
            i = ik;
        end
    end
end

% assignment
X0 = cell(size(X));
foundindex = i>0;
for icol = 1:ncol, X0{icol} = X{icol}(i(foundindex)); end

% final result
if ~isempty(Xfields) % recreate the structure if needed
    X0 = cell2struct(X0,Xfields,1);
elseif reshapeon     % reshape as initial
    X = X0;
    nlig = length(X{1});
    X0 = cell(nlig,ncol);
    for icol = 1:ncol, X0(:,icol) = X{icol}; end
end
if nargout>1, foundindexout = find(foundindex); end