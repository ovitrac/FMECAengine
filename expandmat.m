function [out,outv] = expandmat(object,pattern,dim)
%EXPANDMAT expands an object array along dimension according to a pattern vector (some generalization of repmat)
%   syntax: arrayexpanded = expandmat(array,pattern)
%  options: [arrayexpanded,expansionvector] = expandmat(array,pattern,dim)
%   INPUTS
%           array: numeric, cell, structure array
%         pattern: repetition pattern (array), pattern must be of the same length as the dimension to expand
%             dim: dimension along which repetitions should be applied
%   OUTPUTS
%   arrayexpanded: array with repeated entries along dimension dim
% expansionvector: vector of expansion
%
%   example1: expands 'ABCDE' as 'ABBCCCDDDDEEEEE'
%       txt = 'ABCDE'
%       expandmat(txt,1:length(txt))
%   example2: as example 1 but with a cell array (e contains [1 2 2 3 3 3 4 4 4 4 5 5 5 5 5]')
%       tab = {'A';'B';'C';'D';'E'};
%       [c,e]=expandmat(tab,1:length(tab))
%   example3: duplicates row in a matrix
%       mat = magic(3)
%       m1 = expandmat(mat,1:size(mat,1)) % serialization along dimension 1
%       m2 = expandmat(mat,1:size(mat,2),2) % serialization along dimension 2s

% M2.1 - 07/03/15 - INRA\Olivier Vitrac - rev.

% Revision history

% default
dim_default = 1;

% arg check
if nargin<2, error('two arguments are required'), end
if nargin<3, dim = []; end
if isempty(dim)
    if numel(object)==length(object)
        dim = find(size(object)>1,1,'first');
    end
    if isempty(dim), dim = dim_default; end
end
pattern = pattern(:); npattern = length(pattern);
if dim>ndims(object), error('the dimension (%d) to replicate exceeds the number of dimension %d',dim,ndims(object)), end
if length(pattern)~=size(object,dim)
    error('The size of the pattern (%d) does not match the number of elements along dimension %d',npattern,size(object,dim))
end

% main
patternexpand = arrayfun(@(i,n) i*ones(n,1),(1:npattern)',pattern,'UniformOutput',false);
patternexpand = cat(1,patternexpand{:});
if dim==1,     out = object(patternexpand,:,:,:,:,:,:);
elseif dim==2, out = object(:,patternexpand,:,:,:,:,:);
elseif dim==3, out = object(:,:,patternexpand,:,:,:,:);
elseif dim==4, out = object(:,:,:,patternexpand,:,:,:);
elseif dim==5, out = object(:,:,:,:,patternexpand,:,:);
elseif dim==6, out = object(:,:,:,:,:,patternexpand,:);
elseif dim==7, out = object(:,:,:,:,:,:,patternexpand);
else error('dimension %s is not implemented',dim)
end

% output
if nargout>1, outv = patternexpand; end
