function [X,ndx] = natsortrows(X,col,xpr,varargin)
% Natural-order sort of the rows of a cell array of strings, with customizable numeric format.
%
% (c) 2016 Stephen Cobeldick
%
% Sort the rows of a cell array of strings, sorting the strings by both character
% order and the values of any numeric substrings that occur within the strings.
% The cell of strings must be a matrix (2D). The <col> input of SORTROWS is also
% supported, so NATSORTROWS can be used as a drop-in replacement for SORTROWS.
%
% Syntax:
%  Y = natsortrows(X)
%  Y = natsortrows(X,col)
%  Y = natsortrows(X,col,xpr)
%  Y = natsortrows(X,col,xpr,<options>)
% [Y,ndx] = natsortrows(X,...)
%
% To sort all of the strings in a cell array use NATSORT (File Exchange 34464).
% To sort filenames or filepaths correctly use NATSORTFILES (File Exchange 47434).
%
% See also NATSORT NATSORTFILES SORTROWS SORT CELLSTR NUM2STR NUM2WORDS NUM2ORDINAL NUM2BIP NUM2SIP REGEXP
%
% ### File Dependency ###
%
% This function requires the function NATSORT (File Exchange 34464). The
% inputs <xpr> and <options> are passed directly to NATSORT: see NATSORT for
% case sensitivity, sort direction, numeric substring matching, and other options.
%
% ### Examples ###
%
% A = {'B','2','X';'A','100','X';'B','10','X';'A','2','Y';'A','20','X'};
% sortrows(A) % wrong numeric order:
%  ans = {
%    'A','100','X'
%    'A',  '2','Y'
%    'A', '20','X'
%    'B', '10','X'
%    'B',  '2','X'}
% natsortrows(A) % correct numeric order:
%  ans = {
%    'A',  '2','Y'
%    'A', '20','X'
%    'A','100','X'
%    'B',  '2','X'
%    'B', '10','X'}
% natsortrows(A,'descend')
%  ans = {
%    'B', '10','X'
%    'B',  '2','X'
%    'A','100','X'
%    'A', '20','X'
%    'A',  '2','Y'}
% % Sort ascending by the second column, descending by the third column:
% sortrows(A,[2,-3]) % wrong numeric order:
%  ans = {
%    'B', '10','X'
%    'A','100','X'
%    'A',  '2','Y'
%    'B',  '2','X'
%    'A', '20','X'}
% natsortrows(A,[2,-3]) % correct numeric order:
%  ans = {
%    'A',  '2','Y'
%    'B',  '2','X'
%    'B', '10','X'
%    'A', '20','X'
%    'A','100','X'}
%
% B = {'a','-12';'a','ABCD';'a','3e45';'a','67.8';'a','9';'a','+Inf';'a','NaN'};
% sortrows(B)
%  ans = {
%    'a','+Inf'
%    'a','-12'
%    'a','3e45'
%    'a','67.8'
%    'a','9'
%    'a','ABCD'
%    'a','NaN'}
% natsortrows(B,'ascend','NaN|(+|-)?(Inf|\d+(\.\d+)?((e|E)(+|-)?\d+)?)')
%  ans = {
%    'a','-12'
%    'a','9'
%    'a','67.8'
%    'a','3e45'
%    'a','+Inf'
%    'a','NaN'
%    'a','ABCD'}
%
% C = {'A',10;'a','XX';'A',2;'a',1;'A','#'}; % mixed numeric and string data
% D = cellfun(@num2str,C,'UniformOutput',false); % convert all to string
% natsortrows(D)
%  ans = {
%    'A','#'
%    'a','1'
%    'A','2'
%    'A','10'
%    'a','XX'}
% natsortrows(D,'ascend','\d+','matchcase','asdigit') % see NATSORT for options list
%  ans = {
%    'A','#'
%    'A','2'
%    'A','10'
%    'a','1'
%    'a','XX'}
% natsortrows(D,'ascend','\d+','beforechar') % see NATSORT for options list
%  ans = {
%    'a','1'
%    'A','2'
%    'A','10'
%    'A','#'
%    'a','XX'}
%
% ### Input and Output Arguments ###
%
% Please see NATSORT for a full description of <xpr> and the <options>.
%
% Inputs (*==default):
%  X   = Cell of Strings, with rows to be sorted. Must be 2D (matrix).
%  col = Numeric Vector, column indices to sort <X> by the corresponding
%        columns, where >0=ascending, <0=descending, exactly as per SORTROWS.
%      = String Token, 'descend'/'ascend'* sort direction selection.
%  xpr = String Token, regular expression to detect numeric substrings.
%  <options> can be supplied in any order and are passed directly to NATSORT.
%            Excludes 'descend'/'ascend', which can be supplied as input <col>.
%
% Outputs:
%  Y   = Cell of Strings, input <X> with the rows sorted as per <col>.
%  ndx = Numeric Vector, of size M*1. Row indices such that Y = X(ndx,:).
%
% [Y,ndx] = natsortrows(X,*col,*xpr,<options>)

% ### Input Wrangling ###
%
assert(iscell(X),'First input <X> must be a cell array.')
tmp = cellfun('isclass',X,'char') & cellfun('size',X,1)<2 & cellfun('ndims',X)<3;
assert(all(tmp(:)),'First input <X> must be a cell array of strings (1xN characters).')
assert(ismatrix(X),'First input <X> must be a matrix (size R*C).')
%
if nargin<3
	xpr = '\d+';
end
%
% ### Sort Rows ###
%
[m,n] = size(X);
ndx = (1:m)';
drn = {'descend','ascend'};
%
if nargin<2
	% Sort columns from right to left.
	for k = n:-1:1
		[~,ind] = natsort(X(ndx,k));
		ndx = ndx(ind);
	end
elseif isnumeric(col)
	% Sort columns according to the provided indices.
	assert(isreal(col)&&isvector(col),'Second input <col> must be a real numeric vector.')
	assert(all(floor(col)==col)&&all(abs(col)<=n)&&all(col~=0),...
		'Second input <col> must be a vector of column indices into the first input <X>.')
	varargin{end+1} = [];
	for k = reshape(col(end:-1:1),1,[])
		varargin(end) = drn((3+sign(k))/2);
		[~,ind] = natsort(X(ndx,abs(k)),xpr,varargin{:});
		ndx = ndx(ind);
	end
elseif ischar(col)&&isrow(col)&&any(strcmp(col,drn))
	% Sort all columns descending/ascending.
	varargin{end+1} = col;
	for k = n:-1:1
		[~,ind] = natsort(X(ndx,k),xpr,varargin{:});
		ndx = ndx(ind);
	end
else
	error('Second input <col> must be a numeric vector of indices, or ''ascend''/''descend''')
end
%
X = X(ndx,:);
%
end
%----------------------------------------------------------------------END:natsortrows