function [valid,casout] = checkCAS(cas)
%CHECKCAS return true if a CAS number is valid
% syntax: isvalid = checkCAS(cas) with cas a string or cell array of strings
%         [isvalid,validcas]=checkCAS(cas)
% e.g. [iscas,validcas]=checkCAS('10016-2-3 (alpha) 7585-39-9 (beta)')
%
% For multiple CAS numbers as string:
%   checkCAS('0068515-48-0;0028553-12-0') returns [1 1]
%   checkCAS({'0068515-48-0;0028553-12-0'}) returns 1
%
% See method to check CAS numbers in:
%   http://www.dotnet-news.com/nl/lien.aspx?ID=22056
%   http://www.cas.org/expertise/cascontent/registry/checkdig.html (there is an error)


% MS-MATLAB-WEB 1.0 - 01/05/06 - Olivier Vitrac - rev. 13/07/11

% Revision history
% checkCAS accepts to check multiple CAS numbers as a single string such as checkCAS('0068515-48-0;0028553-12-0')

if iscell(cas)
    valid = false(size(cas));
    for i=1:numel(cas), valid(i) = all(checkCAS(cas{i})); end
elseif ischar(cas)
    valid = false;
    cas = regexp(cas,'([1-9]{1}[0-9]{1,6}-\d{2}-\d)','match');
    if iscell(cas) && length(cas)>1
        valid = checkCAS(cas);
    else
        if ~isempty(cas)
            digits = str2double(regexp(cas{1},'[^-]','match'));
            valid = rem(sum(digits(1:end-1).*(length(digits)-1:-1:1)),10)==digits(end);
        end
    end
else
    valid = false;
end
if nargout>1
    if any(valid), [casout,u] = unique(cas(valid)); valid = valid(u); else casout=''; end
    %if length(casout)==1, casout=casout{1}; end
end
        