function [param,remaining]=argcheck(list,propdefault,keywordlist,varargin)
%ARGCHECK check arguments passed as keywords or as pairs property/value
%    Syntax: param = argcheck(list,propertylist [,keywordlist,'case','property'])
%   Options: [param,remaining] = argcheck(...)
%    INPUTS
%           list: list of arguments coded in a cell array (Note that structure arguments are exapanded in list)
%    propdefault: list of pair property/defaultvalue coded in a structure (propdefault.property=value)
%          propdefault.keyword = defaultvalue
%    keywordlist: list (cell array of strings) of accepted keywords (use [] to omit it)
%         'case': force case sensitive
%     'property': force property precedence (instead of keyword)
%        'order': order fields
%         'keep': keep undefined properties setup in list (e.g. argcheck({'d','D','e','E','f','F'},struct('a','A','b','B','c','C'),'','keep'))
%'nostructexpand': prevent structure to be expanded as a list of parameters
%'NaNequalmissing': NaN are interpreted as missing values
%    OUTPUTS
%          param: recognized keywords and pair property/value (empty values are replaced by their default values)
%            param.keyword = true (if keyword was definded) and false otherwise
%            param.property=value (if property exist), use propertylist.defaultvalue if needed (e.g. when value is set to [])
%      remaining: cell array containing un-used/not recognized parameters
%
% NB: keywords have a higher precedence on properties (a missing keyword is always set to false whatever the value assigned in propdefault)

% MS 2.1 - 25/12/09 - INRA/Olivier Vitrac - rev. 08/04/11

% Revision history
% 27/12/09: case insensitive, several fixes
% 04/01/10: fix ambiguities between keywords and properties
% 01/02/10: fix behavior when no keywords were defined
% 05/02/10: fix syntax typing error
% 09/02/10: fix list as a structure and not a cell
% 09/02/10: fix argcheck([],[],'keyword')
% 09/02/10: add casesensitive
% 12/02/10 fix [a,b]=argcheck({'remaining input'},[],{'cmd','trun'})
% 13/02/10 add 'order'
% 14/02/10 add 'inherit','keep'
% 17/02/10 add 'nostructexpand'
% 24/02/10 fix case for propertylist and when option 'keep' is used
% 08/04/11 add 'NaNequalmissing'

% definitions
options_list = {'case' 'property' 'order' 'keep' 'inherit' 'nostructexpand' 'NaNequalmissing'}; % inherit is a private property

% argcheck
if nargin<1,  error('SYNTAX: [param,remaining]=argcheck(list,propdefault [,keywordlist,''case'',''property''])'), end
if nargin<2, propdefault = []; end
if nargin<3, keywordlist = []; end
if isempty(varargin)
    options = struct('case',false,'property',false,'order',false,'keep',false,'inherit',false,'nostructexpand',false,'NaNequalmissing',false);
else
    options = argcheck(varargin,[],options_list);
end
list = list(:); nlist = length(list);
if isempty(propdefault)
    %provisional implementation (OV 12/02/10)
    %dispf('\tto be used only in DEBUG mode')
    %if (mod(nlist,2)>0), error('list of property/value is unpaired'), end
    if ~isempty(list) && (mod(nlist,2)==0) && iscellstr(list(1:2:end)) 
        if ~iscellstr(list(1:2:end-1)), error('improper key values in list'), end
        propdefault = cell2struct(repmat({[]},nlist/2,1),list(1:2:end-1),1);
    end
end
if options.inherit, param = list; else param = struct([]); end
if ~iscell(list), list = {list}; end

% expand all structures
if ~options.nostructexpand
    istruct = cellfun(@(x)isstruct(x),list);
    nistruct = ~istruct; %find(~istruct);
    istruct = find(istruct);
    nstruct = length(istruct);
    if nstruct>0
        [keys,values] = deal(cell(nstruct,1));
        for i=1:nstruct
            keys{i} = fieldnames(list{istruct(i)});
            values{i} = struct2cell(list{istruct(i)});
        end
        tab = [ cat(1,keys{:}) cat(1,values{:}) ]';
        list = cat(1,list(nistruct),tab(:));
        nlist = length(list);
    end
end

% character inputs
ichar   = find(cellfun(@(x)ischar(x),list));
if ~options.case, list(ichar) = lower(list(ichar)); end
iremaining = true(1,nlist);
    
% search keywords
if ~isempty(keywordlist)
    if ~iscell(keywordlist), keywordlist={keywordlist}; end
    if ~iscellstr(keywordlist), error('list of keywords must be a cell array of char'), end
    for keyword = keywordlist(:)'
        if options.case
            ifound = ismember(list(ichar),keyword{1});
        else
            ifound = ismember(list(ichar),lower(keyword{1}));
        end
        if any(ifound)
            param(1).(keyword{1})=true;
            iremaining(ichar(ifound)) = false;
        else
            param(1).(keyword{1})=false;
        end
    end
end

% search pairs property/value
if ~isempty(propdefault)
    propertylist = fieldnames(propdefault)';
    if options.case, proptosearch = propertylist; else proptosearch = lower(propertylist); end
    if ~iscellstr(propertylist), error('list of properties must be a cell array of char'), end
    iprop = 0;
    for prop = proptosearch
        iprop = iprop + 1;
        ifound = find(ismember(list(ichar),prop));
        if any(ifound)
            for j=1:length(ifound)
                if ~isfield(param,propertylist{iprop}) || options.property
                    if ichar(ifound(j))<nlist && ~(options.NaNequalmissing && any(isnan(list{ichar(ifound(j))+1})))
                        param(1).(propertylist{iprop}) = list{ichar(ifound(j))+1};
                        iremaining(ichar(ifound(j))+[0 1]) = false;
                    elseif j<2
                        param(1).(propertylist{iprop}) =[];
                        iremaining(ichar(ifound(j))) = false;
                    end
                end
            end
        else
            param(1).(propertylist{iprop}) = [];
        end
        if isempty(param(1).(propertylist{iprop})), param(1).(propertylist{iprop}) = propdefault.(propertylist{iprop}); end
    end
end

% keep properties if required
if options.keep
    if options.case % create a structure from remaining parameters
        paramroot = argcheck(list(iremaining),[],'','nostructexpand','case'); % case is used
    else
        paramroot = argcheck(list(iremaining),[],'','nostructexpand'); % no case
    end
    param = argcheck(param,paramroot,'','case','inherit');
end
    
% ouput
if options.order, param = orderfields(param); end
if nargout>1, remaining = list(iremaining); end