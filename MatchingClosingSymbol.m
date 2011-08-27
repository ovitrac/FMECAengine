function pos = MatchingClosingSymbol(str,opensymbol,closingsymbol)
%MATCHINGCLOSINGSYMBOL returns a positions array of matching closing symbols
%   syntax: pos = MatchingClosingSymbol(str [,opensymbol,closingsymbol])
%           str: string
%   opensymbol: default = '(' (any expression is valid, parentheses, braces, curly braces, combination of letters...)
%closingsymbol: default = ')' (any expression is valid)
%           pos: [position of open symbols
%
%   NB: no nesting can be detected when open and closing symbols are identical


% Migration 2.0 - 07/05/2011 - INRA\Olivier Vitrac - rev.


% arg check
if nargin<1, error('one argument is required'); end
if ~ischar(str), error('a string argument is required'); end
if nargin<2, opensymbol = '('; end, opensymbol = regexptranslate('escape',opensymbol);
if nargin<3, closingsymbol = ')'; end, closingsymbol = regexptranslate('escape',closingsymbol);
nonesting =  strcmp(opensymbol,closingsymbol);
    
% search expressions
pos1 = regexp(str,opensymbol)'; 
pos2 = regexp(str,closingsymbol)'; 
if nonesting
    pos1 = pos1(1:2:end);
    pos2 = pos1(2:2:end);
end

% correct sizes
n1=length(pos1);
n2 = length(pos2);
pos1(end+1:end+n2-n1)=NaN; %% add NaN if length(n1)~=length(n2)
pos2(end+1:end+n1-n2)=NaN;
if nonesting, pos = [pos1 pos2]; return; end

% pairing
n1 = length(pos1);
closed = false(n1,1);
order1  = zeros(n1,1);
for i=1:n1
    order1(i) = find((pos2(i)>pos1) & ~closed,1,'last');
    closed(order1(i)) = true;
end
order2(order1) = 1:n1;

% output
pos = [pos1 pos2(order2)];