function [hfindout,shout,repout]=replaceprop(handle,searchrequest,replacerequest)
%REPLACEPROP set properties in graphic objects matching a 
%   SYNTAX replaceprop(handle ,searchrequest,replacerequest)
%          searchrequest: structure or cell array defining the properties to search (all Matlab properties are eligible)
%         replacerequest: structure or cell array setting the new property/values
%       Pair property/value shoulbe either {
%          replaceprop(handle ,searchrequest,struct('property1','value1','property2','value2',...))
%           [hf,search,replace] = replaceprop(...)
%   INPUTS
%       handle: valid handle or array of handles (if empty, search starts from root 0)
%       searchrequest: cell array or structure defining the properties to look for
%       'property1','value1','property2','value2',... replacement properties
%   OUTPUTS
%

% INRA\MS 2.1 - 01/01/2013 - Olivier Vitrac - rev.

% Revision history

% arg check
if nargin~=3, error('3 arguments are required'), end
if ~isempty(handle) && any(~ishandle(handle)), error('some handles are invalid'); end

% search requuest
if isstruct(searchrequest)
    searchrequest = [fieldnames(searchrequest)';struct2cell(searchrequest)'];
    searchrequest = searchrequest(:);
elseif ~iscell(searchrequest)
    error('searchrequest must a structure or a cell array')
end

if rem(length(searchrequest),2), error('the number of pair property/value in searchrequest must be even'), end
n = length(searchrequest)/2;
searchrequest_raw = reshape(searchrequest,2,n);
searchrequest = [searchrequest_raw;repmat({'-and'},1,n)];
searchrequest = searchrequest(:);
searchrequest = searchrequest(1:end-1);

% replacement request
if iscell(replacerequest)
    if rem(length(replacerequest),2), error('the number of pair property/value in searchrequest must be even'), end
    replacerequest = reshape(replacerequest,2,n);
    replacerequest = cell2struct(replacerequest(2,:),replacerequest(1,:),2);
elseif ~isstruct(replacerequest)
    error('replacerequest must a structure or a cell array')
end

% find objects
if ~isempty(handle)
    hfind = findobj(handle,searchrequest{:});
else
    hfind = findobj(searchrequest{:});
end
dispf('REPLACEPROP %d object(s) found',length(hfind))

% do replacement
if any(hfind), set(hfind,replacerequest), end

% outputs
if nargout>0, hfindout = hfind; end
if nargout>1, shout = cell2struct(searchrequest_raw(2,:),searchrequest_raw(1,:),2); end
if nargout>2, repout = replacerequest; end