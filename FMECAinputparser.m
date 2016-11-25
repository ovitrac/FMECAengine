function output = FMECAinputparser(x,varargin)
%FMECAINPUTPARSER generic parser for FMECA inputs
%
% CONSTRUCTOR TO GET THE DEFAULT RULE
%   defaultrule = FMECAinputparser();
%
% SYNTAX TO PARSE A SCALAR INPUT: x, where string is a string or a number
%   xp=FMECAinputparser(x,'property1',value1,'property2',value2,...)
%
% SYNTAX FOR INPUTS VECTORIZED AS STRUCTURES
%   xp=FMECAinputparser(inputs,'rules',rules,'defaultrule',defaultrule)
%   xp=FMECAinputparser(inputs,'rules',rules,'defaultrule',defaultrule'property1',value1,...))
%       with inputs.input1 = x
%            inputs.input1 = [x1 x2 x3...]';
%            inputs.input1 = {x1 x2 x3...}';
%       with 
%            rules.rule1 = struct('property1',value1,'property2',value2,...)
%            rules.rule2 = ...
%
%       Pairing between inputs and rules is done via the property regularexpression
% 
% LIST OF PAIR PROPERTY/VALUE
%                      name: ''
%         regularexpression: ''
%          shortdescription: 'default input'
%           longdescription: 'no description'
%              defaultvalue: []
%                      type: 'numeric'
%                 isliteral: 1
%                       min: -Inf
%                       max: Inf
%                 enablenan: 0
%                 enableinf: 1
%               enableempty: 1
%                ispositive: 1
%                   useunit: 1
%                  quantity: ''
%                      unit: ''
%                    parser: []
%                 errorcode: -100
%              errorsubcode: [1x1 struct]
%              errormessage: ''
%                    islist: 1
%             listseparator: ','
%
% OUTPUTS: output = structure wit fields:
%   value
%   errorcode
%   suberrorcode
%   errormessage
%
%
% EXAMPLE
%{
    defaultrule = struct('type','numeric','ispositive',true,'enablenan',true,'isliteral',false,'useunit',true,'islist',false)
    rules = struct(...
       'length',struct('name','thickness','regularexpression','^l\d*$','shortdescription','l in µm','longdescription','thickness','quantity','l','unit','um','min',.1,'max',1e4), ...
        'Drule',struct('name','D','regularexpression','^D\d*$','shortdescription','D in SI','longdescription','diffusion coefficient','quantity','diff','unit','m2/s','min',1e-22,'max',1e-4),...
         'time',struct('name','t','regularexpression','^t$','shortdescription','t in days','longdescription','time','quantity','t','unit','day','max',10000),...
       'solute',struct('type','text','name','CAS','regularexpression','solute','shortdescription','migrant','longdescription','migrant as defined') ...
        );
    inputs = struct(...
        'l1',{{'1cm' '10' '2304 nm'}},... note that 10 will be interpreted as 10 µm
        'l2',{{'2cm' '20' '1 mm' '1e-4 m'}},...
        'D',{{'4e-12 cm2/s','22e-13'}},...
        't',10 ...
        )
    x = FMECAinputparser(inputs,'defaultrule',defaultrule,'rules',rules)
%}
%
% Migration 2.0 - 31/03/2016 - INRA\Olivier Vitrac - rev. 03/04/2016

% Revision history
% 01/04/2016 first scalar implementation
% 03/04/2016 RC with example

% default
defaultrule = struct(...
    'name','',...
    'regularexpression','',...
    'shortdescription','default input',...
    'longdescription','no description',...
    'defaultvalue',[],...
    'type','numeric',... 'numeric','text','vector','list' 'cell'
    'isliteral',true,...
    'min',-Inf,...
    'max',+inf,...
    'enablenan',false,...
    'enableinf',true,...
    'enableempty',true,...
    'ispositive',true,...
    'useunit',true,...
    'quantity','',...
    'unit','',...
    'parser',[],...
    'errorcode',-100,...
    'errorsubcode',struct('empty',1,'nan',2,'min',3,'max',4,'inf',5,'negative',6,'default',10),...
    'errormessage','',...
    'islist',true,...
    'listseparator',',' ...
    );
default = struct('rules',[],'defaultrule',[]);

%arg check
if nargin<1, output=defaultrule; return, end
[options,unusedpotions] = argcheck(varargin,default,'','nostructexpand');
if ~isempty(unusedpotions)
    options.defaultrule = argcheck(unusedpotions,options.defaultrule,'','keep');
end
if isempty(options.rules), options.rules = struct([]); end
options.defaultrule = argcheck(options.defaultrule,defaultrule);


% vectorization
if isstruct(x)
    inputnames = fieldnames(x); ninputs = length(inputnames);
    rulenames = fieldnames(options.rules); nrules = length(rulenames);
    % pairing inputs and rules
    rules = zeros(ninputs,1);
    if nrules>0
        for irule=1:nrules
             % fix regular expression (if missing, name is taken)
            if ~isfield(options.rules.(rulenames{irule}),'name') || isempty(options.rules.(rulenames{irule}).name)
                options.rules.(rulenames{irule}).name = rulenames{irule};
            end
            if ~isfield(options.rules.(rulenames{irule}),'regularexpression') || isempty(options.rules.(rulenames{irule}).regularexpression)
                options.rules.(rulenames{irule}).regularexpression = sprintf('^%s$',options.strtrim(options.rules.(rulenames{irule}).name));
            end
            % find input(s) matching each prescribed rule (some rules may be also not used)
            inputsfound = ~cellfun(@isempty,regexp(inputnames,options.rules.(rulenames{irule}).regularexpression));
            ialreadyset = find(rules(inputsfound)>0);
            if isempty(ialreadyset)
                rules(inputsfound)=irule;
            else
                dispf('\n%s:: list of conflicting rules')
                dispf('\t%s\n',rulenames{rules([ialreadyset;irule])})
                error('the rules listed above are conflicting, check them and restart %s(...)',mfilename)
            end
        end
    end
    % do [pseudo]recursion
    output = struct([]);
    for iinput=1:ninputs
        if rules(iinput)
            currentrule = {'defaultrule',argcheck(options.rules.(rulenames{rules(iinput)}),options.defaultrule),'rules',''};
        else
            currentrule = {'defaultrule',options.defaultrule,'rules',''};
        end
        if ~iscell(x.(inputnames{iinput})) % scalar input
            output(1).(inputnames{iinput}) = FMECAinputparser(x.(inputnames{iinput}),currentrule{:});
        else % cell inputs
            nx = numel(x.(inputnames{iinput}));
            output(1).(inputnames{iinput}) = repmat(FMECAinputparser(x.(inputnames{iinput}){1},currentrule{:}),nx,1);
            for i=2:nx
                output(1).(inputnames{iinput})(i) = FMECAinputparser(x.(inputnames{iinput}){i},currentrule{:});
            end
        end
    end
    return
elseif ~isempty(options.rules)
    warning('only the defaultrule is used for scalar inputs (rules won''t be used)')
end

% scalar arguments
o = argcheck(options.defaultrule,defaultrule);
description = sprintf('%s (%s)',o.shortdescription,o.name);
defmsg = @(errtxt) [ sprintf('%s has been assigned to the default value: ',description)  errtxt ];
errmsg = @(errtxt,errcode) sprintf('%s cannot be %s (%d)',description,errtxt,errcode);
validparser = isa(o.parser,'function_handle');

% empty
if isempty(x)
    % empty coversion
    switch o.type
        case {'numeric' 'vector'}, xp = [];
        case {'text' 'list'}, xp = '';
        case 'cell', xp={};
    end
    if o.enableempty % empty is accepted
        output = struct('value',xp,'errorcode',0,'suberrorcode',0,'errormessage','');
    elseif any(o.defaultvalue) % empty is not accepted but a default value is proposed
        if ischar(o.defaultvalue) % format the default value
            defaultvaluetxt = o.defaultvalue;
        elseif isnumeric(o.defaultvalue) && length(o.defaultvalue)==1
            defaultvaluetxt = str2double(o.defaultvalue);
        else
            siztxt = arrayfun(@(s) sprintf('%dx',s),size(o.defaultvalue),'UniformOutput',false); siztxt = cat(2,siztxt{:});
            if iscell(o.defaultvalue)
                defaultvaluetxt = sprintf('{%s (%s)}',siztxt(1:end-1),class(o.defaultvalue));
            elseif isstruct(o.defaultvalue)
                defaultvaluetxt = sprintf('[%s struct]',siztxt(1:end-1));
            else
                defaultvaluetxt = sprintf('[%s (%s)]',siztxt(1:end-1),class(o.defaultvalue));
            end
        end
        output = struct('value',o.defaultvalue,'errorcode',0,'suberrorcode',o.errorsubcode.default,'errormessage',defmsg(defaultvaluetxt));
    else % empty is not accepted
        output = struct('value',{{}},'errorcode',0,'suberrorcode',o.errorsubcode.empty,'errormessage',errmsg('empty',o.errrorsubcode.empty));
    end
    return
end

% not empty
if isnumeric(x) % --- numeric data
    switch o.type
        case {'numeric','vector'}
            xp = x;
        case 'text'
            xp = num2str(x);
        case 'list'
            xp = arrayfun(@num2str,x,'UniformOutput',false);
            xp = regexprep(sprintf(['%s' o.listseparator],xp{:}),[o.listseparator '$'],'');
        case 'cell'
            xp = arrayfun(@num2str,x,'UniformOutput',false);
        otherwise
            error('unknown type ''%s''',o.type)
    end
elseif ischar(x)  % --- text data (most of cases)
    if o.isliteral % no interpretation, direct conversion
        switch o.type
            case 'numeric'
                xp = str2double(x);
            case 'vector'
                if o.islist
                    xp = str2double(expandtextaslist(x,o.listseparator));
                else
                    xp = str2double(x);
                end
            case {'text','list'}
                xp = x;
            case 'cell'
                xp = expandtextaslist(x,o.listseparator);
        otherwise
            error('unknown type ''%s''',o.type)
        end
    else % non literal, conversion is required
        if o.islist
            tmp = expandtextaslist(x,o.listseparator);
        else
            tmp = {x};
        end
        for itmp = 1:length(tmp)
            if o.useunit
                tmp{itmp} = FMECAunit(o.quantity,tmp{itmp},'',o.unit);
            elseif validparser
                tmp{itmp} = o.parser(tmp{itmp});
            end
        end % next itmp
        switch o.type % convert to proper type
            case {'numeric' 'vector'}
                xp = [tmp{:}];
            case {'text','list'}
                xp = regexprep(sprintf(['%s' o.listseparator],tmp{:}),[o.listseparator '$'],'');
            case 'cell'
                xp = tmp;
            otherwise
                error('unknown type ''%s''',o.type)
        end
    end % if isliteral
end % if isnumeric    

% outputs (with error messages if needed)
if isnumeric(x)
    if ~o.enablenan && any(isnan(xp))
        output = struct('value',NaN,'errorcode',0,'suberrorcode',o.errorsubcode.nan,'errormessage',errmsg('nan',o.errrorsubcode.nan));
    elseif ~o.enableinf && any(isinf(xp))
        output = struct('value',Inf,'errorcode',0,'suberrorcode',o.errorsubcode.inf,'errormessage',errmsg('Inf',o.errrorsubcode.inf));
    elseif o.ispositive && any(xp<0)
        output = struct('value',NaN,'errorcode',0,'suberrorcode',o.errorsubcode.negative,'errormessage',errmsg('negative',o.errrorsubcode.negative));
    elseif any(xp<o.min)
        output = struct('value',xp,'errorcode',0,'suberrorcode',o.errorsubcode.min,'errormessage',errmsg(sprintf('lower than %0.6g',o.min),o.errrorsubcode.min));
    elseif any(xp>o.max)
        output = struct('value',xp,'errorcode',0,'suberrorcode',o.errorsubcode.max,'errormessage',errmsg(sprintf('greater than %0.6g',o.max),o.errrorsubcode.max));
    else
        output = struct('value',xp,'errorcode',0,'suberrorcode',0,'errormessage','');
    end
else
    output = struct('value',xp,'errorcode',0,'suberrorcode',0,'errormessage','');
end