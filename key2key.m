function [val,dbout,keytreeout] = key2key(db,key,interpreterflag,recursionflag,sample,keynode)
%KEY2KEY returns the values related to keyA:tableA::columnA->tableB::columnB:tableC::columnC->tableD::columnD...
%      SCOPE key2key is a multivalued function (in a Mathematical sense) that implements cardinality operations also known as relationships
%           between tables via primary key, super key and foreign keys. Implemented relationships include:
%           - one-to-one relationships
%           - one-to-many relationships  
%           - many-to-many relationships
%           - cascading between keys
%          Operations above are applied via an intuitive syntax, different of SQL for concision.
%          By default, a UNIQUE operator operator is apply to all relationships involving strings or cell array of strings.
%          The request order is respected only for single requests.
%
%   Examples of relationships in connection with FMECAENGINE
%       PP:polymer::name->classadditives
%       PP:polymer::name->classadditives:substance::class->name
%       PP:polymer::name->classadditives:substance::class->M
%
%       Background:
%           http://en.wikipedia.org/wiki/Many-to-many_(data_model)
%           http://en.wikipedia.org/wiki/Junction_table
%           http://en.wikipedia.org/wiki/Foreign_key
%           http://en.wikipedia.org/wiki/Superkey
%           http://en.wikipedia.org/wiki/Associative_Entities
%
%   Basic syntax: val = key2key(db,key)
%
%  -------------------------------------------------------------------------------
%      INPUTS
%  -------------------------------------------------------------------------------
%             db: structure such as db.tableX.columnY = numerical value or string
%                 alternatively a structure coding for an ODS file such as: db=struct('filename',...,'headers',...)
%            key: string or a 1xn/nx1 string cell array coding for
%                   keyA:tableA::columnA1->columnA2
%                   keyA:tableA::columnA1->columnA2:tableB::columnB1->columnB2
%                   keyA:tableA::columnA1->columnA2:tableB::columnB1->columnB2:tableC::columnC1->columnC2
%            key can contain be combined with numbers or can contain Matlab expressions (see examples for details)
%       NOTE 1: about lists embedded in strings
%            db.tableX.columnY and keys can include list of strings.
%            Values are separated by symbols: ';' (semi-colon) or  '|' (vertical bar)
%            As a rule of thumb, '|' must be preferred in keys as they mimic the behavior or operator OR.
%       NOTE 2: about extra spaces
%            Spaces are tolerated within key values but note that they are trimmed around ';' and '|'.
%            Do not use spaces without checking the reliability of the database content.
%       NOTE 3: about regular expressions
%            Regular expressions are partly implemented (see the rules below).
%            A regular expression must start with symbol '\' and end with '\'
%            Symbol '*' must be replaced by '#' (* is reserved for multiplications).
%            Symbols '(' and ')' must be replaced by '{' and '}' respectively
%            Quantifiers {} must be replaced by {{}}.
%            Quantifier ? is available
%            Operators | . \d \D \s \W are implemented
%            Operators [] ^ $ \< \> (?< (?= (?<! and quantifiers + are not implemented !!!
%       NOTE 4: about multiple keys
%            When multiple keys are combined in expressions, they must be enclosed within ()
%       NOTE 5: about mathematical expressions
%           Almost any valid Matlab expressions can be combined with keys.
%           Diringer() has been implemented within a syntax simpler than with Matlab conventional one
%           Dpiringer(polymercode,key:table1::columnA1->columnB1:table2::columnA2->columnB2 ...->M)
%           Dpiringer(polymercode,key:table1::columnA1->columnB1:table2::columnA2->columnB2 ...->M,temperature)
%           
%  -------------------------------------------------------------------------------
%      OUTPUTS
%  -------------------------------------------------------------------------------
%            val: array of values as stored in db when key is a string
%                 a nx3 array as [min(values median(values) maximum(values)]
%
%   Advanced syntaxes: [val,dbout] = key2key(...)
%                         dbout: database (loaded if needed)
%                      [val,dbout,querytree] = key2key(...)
%                         querytree: structure matching the query and to be used with KEY2KEYGRAPH 
%                           with fields node1, node2,...noden, where nodei is a sub-structure with fields:
%                               parent: parent node (string)
%                                  key: subkey (string)
%                               values: result of the subkey (double array or char array)
%                           isterminal: flag (true if the node is terminal)
%                      NB: in case of multiple subqueries (see example 9), bifurcations are correctly depicted
%                      final results are attached to the longest chain of queries
%               
%   Internal syntax: val = key2key(db,key,interpreterflag,recursionflag)
%                   interpreterflag (default = true) executes keys as possible Matlab expressions
%                   recursionflag (default=false) forces key2key to drop the recursion
%
%   Syntax with specific statistical analysis (works only with cell array of keys)
%                    val = key2key(db,keyarray,[],[],percentiles)
%                    default percentiles = 5th, 50th, 100th
%
%
%  -------------------------------------------------------------------------------
%      EXAMPLES: execute lines between {} with F9 (evaluate)
%  -------------------------------------------------------------------------------
% %{
%  %  EXAMPLE 1 (full toy example)
%        key = 'PP:polymer::name->classadditives:substance::class->name';
%       db.polymer = ...
%           struct('name', {{'LLDPE' 'LLDPE' 'LDPE' 'LDPE' 'MDPE' 'MDPE' 'HDPE' 'HDPE' 'PP' 'PP'}'},...
%        'classadditives', {{'antoxidants' 'antiUV' 'antoxidants' 'antiUV' 'antoxidants' 'antiUV' 'antoxidants' 'antiUV' 'antoxidants' 'antiUV' }'}...
%                 );
%       db.substance = ...
%           struct('class', {{'antoxidants' 'antoxidants' 'antoxidants' 'antiUV' 'antiUV' 'antiUV'}'},...
%                   'name', {{'antiox1' 'antiox2' 'antiox3' 'antiUV1' 'antiUV2' 'antiUV3'}'},...
%                   'M',     [  101       102      103       201       202       203 ]'...
%                   );
%        key2key(db,key)
%        [~,~,details] = key2key(db,key); key2keygraph(details)
% %}
% %{ 
% %  EXAMPLE 2: Variation based on example 1 (not very useful but illustrate the jump between multiple tables)
%        key2 = 'PP:polymer::name->classadditives:substance::class->name:substance::name->M';
%        key2key(db,key2)
%        [~,~,details] = key2key(db,key2); key2keygraph(details)
% %}
% %{
% %  EXAMPLE 3: Variation based on example 2 (use of superclass as key)
%     key3 = 'polyolefin:polymerfamily::name->polymer:polymer::name->classadditives:substance::class->name:substance::name->M';
%     db.polymerfamily = struct(...
%             'name',{{'polyolefin' 'polystyrene'}'},...
%             'polymer',{{'LLDPE;LDPE|MDPE;HDPE;PP' 'HIPS;PS'}'} ....
%             );
%     key2key(db,key3)
%     [~,~,details] = key2key(db,key3); key2keygraph(details)
% %}
% %{
% %  EXAMPLE 4: Variation based on example 3 (behavior in presence of missing values)
%     key4 = 'LDPE|undefined|PP:polymer::name->classadditives:substance::class->name:substance::name->M';
%     key2key(db,key4)
%     [~,~,details] = key2key(db,key4); key2keygraph(details)
% %}
% %{
% %  EXAMPLE 5: Variation based on example 4 (alternative values, note that extra spaces can be included arround operator '|')
%     key5 = 'antiUV2 | antiox1 | antiUV3 :substance::name->M';
%     key2key(db,key5)
%     [~,~,details] = key2key(db,key5); key2keygraph(details)
% %}
% %{
% %  EXAMPLES 6: Variation based on example 5 (search all anitUVi with i=0..9)
%     key6a = '\antiUV\d\ :substance::name->M';
%     key6b = '\antiUV\d{{1}}\ :substance::name->M';
%     key6c = '\antiUV\d{{1,1}}\ :substance::name->M';
%     key6d = '\antiUV\d?\ :substance::name->M';
%     key2key(db,key6a)
%     key2key(db,key6b)
%     key2key(db,key6c)
%     key2key(db,key6d)
%     [~,~,details] = key2key(db,key6a); key2keygraph(details)
% %}
% %{
% %  EXAMPLES 7: Variation based on example 6 (to illustrate the difference beween regular expressions)
%     key7a = '\antiUV{3|2}\ :substance::name->M';
%     key7b = '\antiUV3|2\ :substance::name->M';
%     key7c = '\.#PE.#\:polymer::name->classadditives:substance::class->name:substance::name->M'
%     key7d = '\PE\:polymer::name->classadditives:substance::class->name:substance::name->M'
%     key7e = 'PE:polymer::name->classadditives:substance::class->name:substance::name->M'
%     key2key(db,key7a) % interpreted as containing antiUV3 or containing antiUV2
%     key2key(db,key7b) % interpreted as containing antiUV3 or containing 2
%     key2key(db,key7c)
%     key2key(db,key7d)
%     key2key(db,key7e)
%     [~,~,details] = key2key(db,key7a); key2keygraph(details)
% %}
% %{
% %  EXAMPLES 8: with mathematical formula
%     key8a = 'Dpiringer(PP,PP:polymer::name->classadditives:substance::class->M)'
%     key8b = 'min(Dpiringer(PP,PP:polymer::name->classadditives:substance::class->M,40))';
%     key8c = 'Dpiringer(PP,max(PP:polymer::name->classadditives:substance::class->M))'
%     key8d = 'min(Dpiringer(PP,PP:polymer::name->classadditives:substance::class->M,100))';
%     key8e = 'Dpiringer(LDPE,\antiUV{3|2}\:substance::name->M)';
%     key8f = 'Dpiringer(LDPE,[(\antiUV{3|2}\:substance::name->M);(\antiUV1\:substance::name->M)])';
%     key2key(db,key8a)
%     key2key(db,key8b)
%     key2key(db,key8c)
%     key2key(db,key8d)
%     key2key(db,key8e)
%     key2key(db,key8f)
%    [~,~,details] = key2key(db,key8a); key2keygraph(details)
%    [~,~,details] = key2key(db,key8f); key2keygraph(details)
% %}
% %{
% %  EXAMPLE 9: example 8a with double and triple keys
%     db.contact = struct('condition',{{'cond1' 'cond2' 'cond3'}},'temperature',[25 40 110]);
%     key9a = 'Dpiringer(PP,(PP:polymer::name->classadditives:substance::class->M),(cond3:contact::condition->temperature))';
%     key2key(db,key9a)
%     db.scenario = struct('id',{{'scenarioA' 'scenarioB' 'scenarioC'}},'polymer',{{'PP' 'LDPE' 'PS'}});
%     key9b = 'Dpiringer((scenarioA:scenario::id->polymer),(PP:polymer::name->classadditives:substance::class->M),(cond3:contact::condition->temperature))';
%     key2key(db,key9b)
%     key9c = 'Dpiringer((scenarioA:scenario::id->polymer),(scenarioA:scenario::id->polymer:polymer::name->classadditives:substance::class->M),(cond3:contact::condition->temperature))';
%     key2key(db,key9c)
%     [~,~,details] = key2key(db,key9a); key2keygraph(details)
%     [~,~,details] = key2key(db,key9c); key2keygraph(details)
% %}
%
% %  EXAMPLE 10: with list contains number as CAS number
%     db.substance =  ...
%           struct('class', {{'antoxidants' 'antoxidants' 'antoxidants' 'antiUV' 'antiUV' 'antiUV'}'},...
%                   'name', {{'antiox1' 'antiox2' 'antiox3' 'antiUV1' 'antiUV2' 'antiUV3'}'},...
%                   'M',     [  101       102      103       201       202       203 ]',...
%                   'CAS', {{'0000112-84-5' '0002082-79-3' '0006683-19-8' '0041484-35-9' '0000123-28-4' '0000057-11-4'}'},...
%                   'CASm', {{'C0000112-84-5' 'C0002082-79-3' 'C0006683-19-8' 'C0041484-35-9' 'C0000123-28-4' 'C0000057-11-4'}'});
%    db.result = struct('sample', {{'S1' 'S2'}'},...
%                       'mol',{{'0000112-84-5;0000057-11-4;0041484-35-9' '0002082-79-3;0006683-19-8'}'},...
%                       'molm',{{'C0000112-84-5;C0000057-11-4;C0041484-35-9' 'C0002082-79-3;C0006683-19-8'}'});
%     key10a = 'S1:result::sample->mol:substance::CAS->name'; 
%     key2key(db,key10a)
%     key10b = 'S1:result::sample->molm:substance::CASm->name'; 
%     [~,~,details] = key2key(db,key10b); key2keygraph(details)
%
%   See also: FMECAENGINE FMECASINGLE ISMEMBERLIST EXPANDTEXTASLIST LOADFMECAENGINEDB KEY2KEYGRAPH


% Migration 2.0 - 06/05/2011 - INRA\Olivier Vitrac - rev. 29/12/11

% Revision history
% 07/05/11 add Matlab expressions
% 08/05/11 try to interpret literally invalid keys as Matlab expression
% 10/05/11 enable a database filename and literal numbers
% 17/07/11 replace ismember by ismemberlist, update help and expand examples accordingly
% 18/07/11 enable symbols # { } ? and , (except if Dpiringer is detected), improved help
% 25/07/11 add examples 9 (multiple keys within Dpiringer)
% 25/08/11 split final results that include lists of strings (with ';' or '|')
% 26/08/11 fix isDpiringer for cell arrays, add sample
% 23/12/11 add KEYTREE, keytreeout
% 27/12/11 do not update KEYTREE when keytreeout is not requested
% 27/12/11 keytreeout returns val and isterminal
% 28/12/11 modified examples to plot the the query tree
% 29/12/11 fix multiple keys, fix static/dynamic nodes
% 07/01/15 add example 10 with CAS number list for debugging

%default (do not modify the regular expressions without keeping a copy)
% patternkey = '^(?<key>[A-Za-z0-9_\.\-;,\|\s\\#\{\}\?]+):(?<tablesource>[A-Za-z][A-Za-z0-9_]*)::(?<columnsource>[A-Za-z][A-Za-z0-9_]*)->(?<columndestination>[A-Za-z][A-Za-z0-9_]*)';
patternkey = '^(?<key>[A-Za-z\\][A-Za-z0-9_\.\-;,\|\s\\#\{\}\?]*):(?<tablesource>[A-Za-z][A-Za-z0-9_]*)::(?<columnsource>[A-Za-z][A-Za-z0-9_]*)->(?<columndestination>[A-Za-z][A-Za-z0-9_]*)';
% patternisakey = '([a-zA-Z0-9_\.\-;,\|\s\\#\{\}\?]+:.*->[A-Za-z][A-Za-z0-9_]*)';
% patternisakey = '([a-zA-Z0-9_\.\-;,\|\s\\#\{\}\?]+:[^\(\);]*->[A-Za-z][A-Za-z0-9_]*)'; %to enable a separation of successive keys
patternisakey = '([A-Za-z\\][a-zA-Z0-9_\.\-;,\|\s\\#\{\}\?]*:[^\(\);]*->[A-Za-z][A-Za-z0-9_]*)'; % key must start by a character or \
num = '([+-]?\s*\d+\.?\d*[eEdD]?[+-]?\d*)'; %regular expression for any scientific numbers
splitstr = '\s*[;\|]+\s*';
interpreterflag_default = true;
recursionflag_default = false;
sample_default = [5 50 95];
nodename = 'node%d';
rootnode = sprintf(nodename,1);

% Peristent variable
persistent KEYTREE

% arg check
if nargin<2, error('Two arguments are required'); end
if ~isstruct(db), error('the first argument must be a structure'); end
if isempty(key), val = []; if nargout>1, dbout=db; end; if nargout>2, keytreeout=struct([]); end, return; end
if nargin==2, if nargout>2, KEYTREE = struct('lastidx',0); else KEYTREE = struct([]); end, end
if nargin<3, interpreterflag = []; end
if nargin<4, recursionflag = []; end
if nargin<5, sample = []; end
if nargin<6, keynode = []; end
if isempty(interpreterflag), interpreterflag = interpreterflag_default; end
if isempty(recursionflag), recursionflag = recursionflag_default; end
if isempty(sample), sample = sample_default; end
if isempty(keynode), keynode = struct('parent','','nodeidx',0); end
sample = sample(:)';
useKEYTREE = isfield(KEYTREE,'lastidx');

% Remove , when Dpiringer is used
if ischar(key) 
    isDpiringer = any(strfind(key,'Dpiringer'));
elseif iscellstr(key)
    isDpiringer  = any(~cellfun('isempty',regexp(key,'Dpiringer')));
else
    error('unknown key type')
end
if isDpiringer %any(strfind(key,'Dpiringer'))
    patternkey(patternkey==',')='';
    patternisakey(patternisakey==',')='';
end
    
% Recursion
if ~recursionflag
    if isnumeric(key) % already numeric: no conversion
        val = key;
        if nargout>1, dbout=db; end;
        return
    elseif iscell(key) % cell array containing keys
        nkey = length(key);
        tmp = key2key(db,key{1},interpreterflag_default); % interpretation flag set to true (as default) to avoid resetting KEYTREE
        if isnumeric(tmp)
            val =zeros(nkey,3); 
            for i=1:nkey
                if isnumeric(key{i}) % keep the value
                    val(i,:) = [1 1 1]*key{i};
                else % key to interpret
                    tmp = key2key(db,key{i},interpreterflag_default); % interpretation flag set to true (as default) to avoid resetting KEYTREE
                    % before sample
                    %nnan = ~isnan(tmp);
                    %if any(nnan), med = median(tmp(nnan)); else med = NaN; end
                    %val(i,:) = [min(tmp) med max(tmp)]; % keep only min median and max
                    val(i,:) = prctile(tmp(~isnan(tmp)),sample);
                end
            end
        else
            val =cell(nkey); for i=1:nkey, val{i} = key2key(db,key{i},interpreterflag_default); end % interpretation flag set to true (as default) to avoid resetting KEYTREE
        end
        if nargout>1, dbout=db; end;
        if nargout>2
            if isfield(KEYTREE,rootnode), KEYTREE.(rootnode).values = val; end
            keytreeout = addisterminal(KEYTREE);
        end
        return
    elseif ischar(key) && interpreterflag && ~isempty(regexp(key,patternisakey,'once')) % MULTIPLE KEYS
        %%% convert the key expression, add terminal ';', remove assignment if any
        staticnode = true; % static by default (no bifurcation)
        if useKEYTREE
            if keynode.nodeidx>0, parent = sprintf(nodename,keynode.nodeidx); else parent=''; end
            if length(regexp(key,patternisakey))>1 || isDpiringer % multiple keys => dynamic nextkeynode => bifurcatons
                staticnode = false;
                KEYTREE(1).(sprintf(nodename,KEYTREE.lastidx+1)) = struct('parent',parent,'key',key,'values',[]); % store the key in KEYTREE
                KEYTREE.lastidx = KEYTREE.lastidx + 1;
                nextkeynode = @(currentidx) struct('parent',sprintf(nodename,KEYTREE.lastidx),'nodeidx',currentidx+1,'values',[]); %#ok<NASGU>
                key = regexprep(key,{patternisakey,'([^;])$','^[^(]*='},{'key2key(db,''$1'',false,false,sample,nextkeynode(KEYTREE.lastidx))','$1;',''});
            else % single key
                nextkeynode = struct('parent',parent,'nodeidx',keynode.nodeidx+1,'values',[]); %#ok<NASGU>
            end
        else % static nextkeynode
            nextkeynode = struct([]); %#ok<NASGU>
        end        %%% check whether Diringer is used
        if staticnode % default behavior
            key = regexprep(key,{patternisakey,'([^;])$','^[^(]*='},{'key2key(db,''$1'',false,false,sample,nextkeynode)','$1;',''});
        end

        % indentify the span of Dpiringer(...) even if it contains subexpressions with additional ()
        [start,stop] = regexpi(key,'Dpiringer\(');
        if length(stop)>1, error('Error in ''%s''\nMultiple expressions of Dpiringer() are not allowed (no restrictions for other functions).',key), end 
        if ~isempty(start)
            bracepos = MatchingClosingSymbol(key);
            stop = bracepos(bracepos(:,1)==stop,2);
            % protects the first argument (coding for the polymer name)
            subkey = key(start:stop);
            [arg1,~,stop1] = regexp(subkey,'Dpiringer\((.*?),','tokens'); arg1 = uncell(arg1);
            if arg1{1}(1)~='''' && isempty(regexp(arg1{1},'\s*\(?\s*key2key\(', 'once')) % nb: polymer can be given by a key
                arg1{1} = ['''' arg1{1} ''''];
            end
            key = [key(1:start-1) sprintf('Dpiringer(%s%s',arg1{1},subkey(stop1:end)) key(stop+1:end)];
        end %Dpiringer
        % check whether the database must be loaded
        if isfield(db,'filename') && isfield(db,'headers'), db = loadods(db.filename,'sheetname','all','headers',db.headers); end
        try
            val = eval(key);
        catch err
            dispf('Error with the key ''%s''',key)
            rethrow(err)
        end
        if nargout>1, dbout=db; end;
        if nargout>2
            if isfield(KEYTREE,rootnode), KEYTREE.(rootnode).values = val; end
            keytreeout = addisterminal(KEYTREE);
        end
        return
    end
end % recursion flag

%%% execute key code
if useKEYTREE
    KEYTREE(1).(sprintf(nodename,keynode.nodeidx)) = struct('parent',keynode.parent,'key',key,'values',[]); % store the key in KEYTREE
    KEYTREE.lastidx = KEYTREE.lastidx + 1;
end
[start,stop] = regexp(key,patternkey);
if isempty(start) % invalid key but try to interpret literally the code
    try
        val = eval(key);
        if nargout>1, dbout=db; end;
        if useKEYTREE, KEYTREE.(sprintf(nodename,keynode.nodeidx)).val = val; end
        return
    catch %#ok<CTCH>
        dispf('WARNING: unrecognized syntax for ''%s''.\nTry to read numbers literally, please check',key)
        val = cellfun(@(x)str2double(x),regexp(key,num,'match'));
        if isempty(val)
            error('unable to understand the key code ''%s'', it does not return/contain any numbers',key);
        end
        if nargout>1, dbout=db; end;
        return
    end
end

%valid key
remainderkey  = key(stop+1:end);
op            = regexp(key,patternkey,'names');
if ~isfield(db,op.tablesource), error('the source table ''%s::'' is missing in ''%s''',op.tablesource,key); end
if ~isfield(db.(op.tablesource),op.columnsource), error('the source column ''%s::%s'' is missing in ''%s''',op.tablesource,op.columnsource,key); end
if ~isfield(db.(op.tablesource),op.columndestination), error('the destination column ''%s::%s'' is missing in ''%s''',op.tablesource,op.columndestination,key); end

% read columnsource (see ismemberlist for details, by default: keepBorder,regularexpr are set to TRUE)
%found = db.(op.tablesource).(op.columndestination)(ismember(db.(op.tablesource).(op.columnsource),op.key));
found = db.(op.tablesource).(op.columndestination)(ismemberlist(db.(op.tablesource).(op.columnsource),op.key,splitstr,true,true));
if iscellstr(found)
    found = unique(found);
elseif iscell(found)
    if all(cellfun(@(x) isnumeric(x), found))
        found = cat(1,found{:});
    end
end

% outputs
if isempty(remainderkey)
    
    if iscellstr(found) % split results that include list
        val = uncell(regexp(found,splitstr,'split'));
        if ~iscellstr(val)
            for i=1:length(val), if size(val{i},2)>1, val{i} = val{i}'; end, end
            val = cat(1,val{:});
        else
            if size(val,2)>1, val = val'; end
        end
    else
        val = found;
    end
    if useKEYTREE, KEYTREE.(sprintf(nodename,keynode.nodeidx)).values = val; end
        
else
    
    if useKEYTREE
        nextkeynode = struct('parent',sprintf(nodename,keynode.nodeidx),'nodeidx',[]);
    else
        nextkeynode = struct([]);
    end
    val = cell(0,1);
    for i=1:length(found)
        if useKEYTREE, nextkeynode.nodeidx = KEYTREE.lastidx+1; end
        tmp = key2key(db,sprintf('%s%s',found{i},remainderkey),false,true,sample,nextkeynode); % call with recursion
        if i==1, val = tmp; else val = [val;tmp]; end %#ok<AGROW>
        if useKEYTREE
            KEYTREE.(sprintf(nodename,nextkeynode.nodeidx)).values = tmp;
        end
    end
    
end

%% additional outputs
if nargout>1, dbout = db; end
if nargout>2
    if isfield(KEYTREE,rootnode), KEYTREE.(rootnode).values = val; end
    keytreeout = addisterminal(KEYTREE);
end

end % end main function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIVATE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDISTERMINAL set isterminal field to true when the node is terminal
function keytree = addisterminal(keytree)
    keytree = rmfield(keytree,'lastidx');
    nodetreelist = fieldnames(keytree)';
    paths = buildmarkov(keytree);
    terminalnodes = cellfun(@(f) f{end},paths,'UniformOutput',false); %unique not required as (single parent)
    for f = nodetreelist
        keytree.(f{1}).isterminal = ismember(f,terminalnodes);
        % force all nodes with numeric data to be terminal (GOAL: forces key2keygraph to display the results of intermediate queries)
        if isnumeric(keytree.(f{1}).values) && ~isempty(keytree.(f{1}).values)
            keytree.(f{1}).isterminal = true;
        end
    end
end % end keytree