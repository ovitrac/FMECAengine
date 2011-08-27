function [dupval,dupind,ndup] = findduplicates(x,forcenanequal)
%FINDDUPLICATES finds duplicate values in vectors (numeric array or cell array of strings)
%   [dupval,dupind,ndup] = findduplicates(x [,forcenanequals])
%           dupval: duplicated values
%           dupind: indices of duplicated values
%          ndup: number of repetitions
% forcenanequal: flag to force NaN as a single value (default=true)
%
%   See also: find_multiple

% MS 2.1 - 12/04/2011 INRA\Olivier Vitrac - rev. 13/05/2011

% Revision history
% 13/05/2011 force 0x0 cell (instead of 0x1 cell)

% arg check
if nargin<1, error('one argument is required'); end
if nargin<2, forcenanequal = []; end
if ~isnumeric(x) && ~iscellstr(x), error('x must be a numeric array or cell array of string'); end
if isempty(forcenanequal), forcenanequal = true; end
forcenanequal = forcenanequal && isnumeric(x);
   

% compute
[xu,iu,ju] = unique(x(:)); %#ok<ASGLU>
if forcenanequal
    iNaN = isnan(xu);
    if any(iNaN)
        ifirstNaN = find(iNaN,1,'first');
        iNaN(ifirstNaN) = false; % keep first only
        xu = xu(~iNaN);
        ju(ju>length(xu)) = ifirstNaN;
    end
end

found = find_multiple(1:length(xu),ju,1); % counts xu values
counts = cellfun(@(x) length(find(x>0)),num2cell(found,1))'; % alternative: counts = histc(ju,1:length(xu))
dupval = xu(counts>1);
if isempty(dupval), dupval = {}; end

% additional outputs
if nargout>1, idx = (1:numel(x))'; dupind = idx(counts(ju)>1); end
if nargout>2, ndup = counts(counts>1); end