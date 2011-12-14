function deleteplotpub(hplot)
% DELETEPLOTPUB delete PLOTPUB objects
%   Syntax: deleteplotpubs(hplot)
%
%   see also: PLOTPUB

% MS 2.0 - 16/08/07 - Olivier Vitrac - rev. 31/08/07

% history
% 18/08/07 fix empty handles
% 30/08/07 fix delete text objects
% 31/08/07 replace length by numel

% arg check
if nargin~=1, error('one argument is required'), end
if ~isstruct(hplot) || ~isfield(hplot,'leg') || ~isfield(hplot,'line') || ~isfield(hplot,'marker') || ~isfield(hplot,'text')
    error('invalid handle object, hplot must be created with PLOTPUB')
end

% delete all handles
m = numel(hplot);
for f=fieldnames(hplot)'
    for i=1:m
        h = hplot(i).(f{1});
        if any(h) && all(ishandle(h)), delete(h); end
    end
end