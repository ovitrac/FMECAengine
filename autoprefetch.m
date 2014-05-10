function out = autoprefetch(varargin)
%AUTOPREFETCH generates and manages automatically PREFETCH files
%
%NOTE1: self-contained OBJECT oriented code, implemented as Javascript does for robustness
%NOTE2: The main intent is the acceleration of any piece of Matlab code by reusing one or several prefetch files.
%NOTE3: 3 versions of the same prefetech file are automatically managed by ROTATIONMAT
%NOTE4: Several sections of a same same code can be managed with different prefetch files (mblock is incremented).
%NOTE5: The database (managed via persistent variables) is shared between Matlab codes attached to the same Matlab Session
%NOTE6: A cache of variables is managed in memory, use memory to check the amount of physical memory allocated to Matlab.
%       Clear all will free memory.
%
%  CONTRUCTORS/SYNTAX p=autoprefetch('keyword',...)
%   p=autoprefetch('open'...)
%   p=autoprefetch('append'...)
%   'open': reset the internal database (in memory) and add a new session
%   'append': add a new session to the existing database
%
%  LIST OF IMPLEMENTED METHODS/SYNTAX p.method(...) or out=p.method(...)
%  USE TAB after p. to list all available methods
%       p.read(...) lists the variables currently used
%                   returns TRUE if the file exist (to be used within IF) or if 'forcedrefresh' is used
%     p.update(...) calculates the diff between now and the variables listed with read()
%                   returns refreshed variable names (to be used after read())
%      p.write(...) saves the variables (in memory and in a prefetch file) refreshed with update()
%                   returns the name of the prefetch file
%                   use 'overwrite' to force overwriting on existing prefetch or global files
%      p.close(...) removes the curent session from memory
%                   returns true
%       p.save(...) save all sessions in a same global gile
%                   returns the name of the file
%   p.closeall(...) closes all the internal databases
%                   returns true
%       p.disp(...) displays statistics on the current session
%                   returns the object definition
%       p.help(...) displays the available attributes to be used as (...)
%                   the same attributes can be set at construction time
%                   NOTE: all set attributes modify permamently the created object
%
%     --- LIST OF PUBLIC PROPERTIES/ATTRIBUTES ---
%     [   'keyword'] [type   ] [description                                  ] [default value   ]
%     [------------] [-------] [---------------------------------------------] [----------------]
%     [   mfilename] [string ] [name of the M file                           ] [mfilename()     ]
%     [      mblock] [integer] [block index of the Mfile                     ] [last block +1   ]
%     [        path] [string ] [default path                                 ] [''              ]
%     [prefetchfile] [string ] [name of the prefetch file                    ] [AUTOPREFETCH_   ]
%     [  globalfile] [string ] [name of the global file                      ] [AUTOPREFETCH_ALL]
%     [   frequency] [scalar ] [frequency of file rotation in days           ] [1               ]
%     [    nrotates] [integer] [number of copies of prefetch and global files] [2               ]
%     [   rotateeco] [integer] [if true only NULL files are saved at creation] [1               ]
%     [   overwrite] [integer] [if true all files will be overwritten        ] [false           ]
%     --- all keywords can be used for constructors and methods ---
    
%  TYPICAL USAGE IN A SCRIPT/FUNCTION
%     p = autoprefetch('open',....) % creates the session at the beginning of the file
%     [...]
%     if p.read()     % check whether the prefetch file exist or not, load data, if not it loads it
%           [...]
%           p.write()
%     end
%     [...]
%     p = autoprefetch('append',....)
%     if p.read()...
%           [...]
%           p.write()
%     end
%     [...]
%    p.close()
%

% INRA\MS 2.1 - 24/12/12 - Olivier Vitrac - rev. 31/12/12

% Revision history
% 31/12/12 release candidate, with full help
% 01/01/13 fix forcedrefresh and all flags to be inherited, add overwite and forcedrefresh in disp()
% 02/01/12 cosmetic changes in the code, implementation of automatic mlmfile



%% todo: OPEN/CLOSE 

% PERSISTENT variables
% Session data are 
% STATE.session stores all session data (session = one single prefetch)
% MFILES(i).mfilename and MFILES(i).mblock
persistent STATE MFILES
if isempty(MFILES) && ~isstruct(MFILES)
    MFILES = repmat(struct('mfilename','','mblock','','session',''),0,1);
end
if isempty(STATE) && ~isstruct(STATE), STATE = struct([]); end


% default
methodflags = {'open','append','read','update','write','close','closeall','save','disp','help'};
userflags = {'forcedrefresh','overwrite'};
keywords = [methodflags userflags];
default = struct(...
    'mfilename','',...
    'mblock',[],...
    'path',[],...
    'forcedrefresh',false,...
    'prefetchfile','',...
    'globalfile','',...
    'frequency',1,... 1 days
    'nrotates', 2, ...
    'rotateeco',true,...
    'overwrite',false,...
    'open',false,...
    'append',false,...
    'read',false,...
    'update',false,... private method (internally used)
    'write',false,...
    'close',false,...
    'closeall',false,...
    'save',false,...
    'disp',false,...
    'help',false,...
    'session','',...
    'inherit',[],...
    'varnames',{{}},...
    'vars',{{}},...
    'action1','',...
    'action2','',...
    'override',{{}} ...
    );
allproperties = fieldnames(default);
% user properties (to be modfied by user methods)
publicproperties = {'mfilename';'mblock';'path';'prefetchfile';'globalfile';'frequency';'nrotates';'rotateeco';'overwrite';'forcedrefresh'};
% system properties
privateproperties = setdiff(allproperties,publicproperties);
% inheritable properties
inheritableproperties = [publicproperties;{'session'}];

% arg check
o = argcheck(varargin,default,keywords); % user arguments/parameters
if ~isempty(o.inherit), o = argcheck(o.inherit,rmfield(o,'inherit'),'','keep','case'); end %inheritence from object (without methods)

% user override (if any) of public properties
% to be used with user methods: read(), update(), write(), save(),...
sessionok = ~isempty(o.session) && isfield(STATE,o.session);
if ~isempty(o.override) && sessionok
    override = argcheck(o.override,rmfield(o,intersect(fieldnames(o),privateproperties)),userflags);
    if isempty(override.prefetchfile), error('AUTOPRFETCH: empty prefetch filename passed as user override'), end
    override.prefetchfile = validprefetchfile(override.prefetchfile,override.path,override.overwrite || override.forcedrefresh,override.mfilename);
    STATE.(o.session).override = override;
    o.override = {};
end
if sessionok && ~isempty(STATE.(o.session).override)
    o = argcheck(STATE.(o.session).override,o);
end

% ###################################################################################################
%       ACTIONS DEFINED VIA FLAGS
% ###################################################################################################

if o.append || o.open
    % set a valid mfile name
    %if isempty(o.mfilename), o.mfilename = evalin('caller','mfilename'); end
    if isempty(o.mfilename), o.mfilename = 'commandline'; end
    % create a new session
    if isempty(o.session), o.session = autoprojectname('',true); end
    imf = find(ismember({MFILES.mfilename},o.mfilename));
    if o.open %% fully new session
        if isempty(o.mblock), o.mblock = 1; end
        if ~isempty(imf)
            dispf('AUTOPREFETCH: WARNING delete previous sessions (forced close) for file ''%''',o.mfilename)
            STATE = rmfield(STATE,intersect(fieldnames(STATE),{MFILES(imf).session}));
            MFILES(imf) = [];
        end
        dispf('AUTOPREFETCH open a new session for file ''%s'' - block=%d',o.mfilename,o.mblock)
    else % append to previous sessions
        % search for existing mfilename
        if length(imf)<1, error('the M file ''%s'' is not in the database, use ''open'' instead of ''append''',o.mfilename), end
        o.mblock = max([MFILES(imf).mblock])+1; % increment the last existing block
        dispf('AUTOPREFETCH append to an existing session for file ''%s'' - block=%d',o.mfilename,o.mblock)
    end
    MFILES(end+1) = struct('mfilename',o.mfilename,'mblock',o.mblock,'session',o.session);
    % update the prefetchfile
    %if isempty(o.prefetchfile), o.prefetchfile = sprintf('AUTOPREFETCH_%s_block%02d_%s.mat',o.mfilename,o.mblock,datestr(now,'yyyy-mm-dd_HH_MM')); end
    if isempty(o.prefetchfile), o.prefetchfile = sprintf('AUTOPREFETCH_%s_block%02d.mat',o.mfilename,o.mblock); end
    o.prefetchfile = validprefetchfile(o.prefetchfile,o.path,true,o.mfilename); % force overwrite at creation (no error if file exist)
    % add indentifiers to the new session
    nfo = struct(...
        'date',datestr(now,'yyyy-mm-dd HH:MM'),...
        'localname',localname,...
        'path',pwd,...
        'mfilename',o.mfilename,...
        'prefetchfile',o.prefetchfile,...
        'mblock',o.mblock,...
        'session',o.session ...
        );
    STATE(1).(o.session) = struct('varlist',[],'varnames',[],'vars',[],'nfo',nfo,'override',[]);
    
    % ###################################################################################################
    %         CREATE USER METHODS (2 possibilities so-called CASE 1 and CASE 2)
    % ###################################################################################################
    %
    % NOTE: auto cleanup not working due to a bug in Matlab
    % cleanup = sprintf('onCleanup(@() autoprefetch(struct(''close'',true,''override'',{{}},''inherit'',struct(''session'',''%s''))))',o.session);
    %
    % # CASE 1 (the prefetchfile does not exist or forcedrefresh is used)
    if ~exist(o.prefetchfile,'file') || o.forcedrefresh
        read = '{substructarray(whos,''name'')}';
        actionread = {'[]';'[]'}; %cleanup
    else % # CASE 2 (the prefetch file exist, we will use it)
        read = sprintf('{setdiff(substructarray(whos,''name''),substructarray(whos(''-file'',''%s''),''name''))}',o.prefetchfile);
        actionread = {
            sprintf('cprintf(''%%s\\n'',''AUTOPREFETCH: [%s - %d] loading the following prefetch file...'',fileinfo(''%s'',[],false,true))', ...
            o.mfilename,o.mblock,o.prefetchfile); ...  first action = some displays (returned as text, read() will display it)
            sprintf('[evalc(''load(''''%s'''')'') '''']',o.prefetchfile) % second action = loading data with a capture to get a text output
            };
    end
    oinherit = rmfield(o,setdiff(allproperties,inheritableproperties)); % keep only inheritable properties
    
    % # METHOD read()
    % the method read must return TRUE (to be used within IF)
    out.read =  @(varargin) autoprefetch(...
        struct('read',true,... flag of the method (read here means read/collecting the list of variables)
        'varnames',evalin('caller',read),... main action = retrieving existing variables
        'action1',evalin('caller',actionread{1}),... action to be performed
        'action2',evalin('caller',actionread{2}),... action to be performed         'action3',{evalin('caller',actionread{3})},... action to be performed
        'override',{varargin},... user override properties
        'inherit',oinherit) ); % inheritance properties
    
    % # METHOD update()
    %  refresh variable names (to be used after read())
    update = '{substructarray(whos,''name'')}';
    out.update =  @(varargin) autoprefetch(...
        struct('update',true,... flag of the method (read here means read/collecting the list of variables)
        'varnames',evalin('caller',update),... main action = retrieving existing variables
        'override',{varargin},... user override properties
        'inherit',oinherit) ); % inheritance properties
    
    % # METHOD write()
    %    => step1: refresh variable names (call update() first)
    %    => step2: copy all data as a structure
    update = sprintf('struct(''update'',true,''varnames'',{substructarray(whos,''name'')},''inherit'',struct(''session'',''%s'',''mblock'',%d,''mfilename'',''%s'',''prefetchfile'',''%s''))',...
        o.session,o.mblock,o.mfilename,o.prefetchfile);
    out.write =  @(varargin) autoprefetch(...
        struct('write',true,... flag of the method
        'vars', {evalin('caller',sprintf('{%s}''',cprintf(' %s ',evalin('caller',sprintf('autoprefetch(%s)',update)))))},... retrieve the variables content
        'override',{varargin},... user override properties
        'inherit',oinherit) ); % inheritance properties
    
    % # METHOD close()
    out.close =  @(varargin) autoprefetch(...
            struct('close',true,... flag of the method
                   'override',{varargin},... user override properties
                   'inherit',oinherit) ); % inheritance properties
    
    % # METHOD closeall()
    out.closeall =  @(varargin) autoprefetch(...
            struct('closeall',true,... flag of the method
                   'override',{varargin},... user override properties
                   'inherit',oinherit) ); % inheritance properties
               
    % # METHOD save()
    out.save =  @(varargin) autoprefetch(...
            struct('save',true,... flag of the method
                   'override',{varargin},... user override properties
                   'inherit',oinherit) ); % inheritance properties
               
    % # METHOD disp()
    out.disp =  @(varargin) autoprefetch(...
            struct('disp',true,... flag of the method
                   'override',{varargin},... user override properties
                   'inherit',oinherit) ); % inheritance properties
               
    % # METHOD help()
    out.help =  @(varargin) autoprefetch(...
            struct('help',true,... flag of the method
                   'override',{varargin},... user override properties
                   'inherit',oinherit) ); % inheritance properties
               
    % ###################################################################################################
    %         APPLY METHODS
    % ###################################################################################################
    
elseif o.read
    % ### read() method
    % # the method list the variables currently used
    % # the method returns true if the prefetch file does not exist or if forcedrefresh is used
    if ~isempty(o.action1) && ischar(o.action1), disp(o.action1); end
    if ~isempty(o.action2) && ischar(o.action2), disp(o.action2); end
    %if ~isempty(o.action3) && ischar(o.action3), disp(o.action3); end
    validsession('read()')
    if ~isempty(STATE.(o.session).varnames) || ~isempty(STATE.(o.session).vars)
        dispf('\nAUTOPREFETCH read() error: corrupted call for the following prefetch session:\n\t--DUMP --')
        STATE.(o.session).nfo
        dispf('\t--DUMP --\n')
        error('method read() can be used only once (see above), use: p = autoprefetch(''append'',....) to create a new session')
    elseif ~isempty(STATE.(o.session).varlist)
        dispf('AUTOPREFETCH: update the list of variables (not to be saved) for the file ''%s'', block=%d',o.mfilename,o.mblock)
    else
        dispf('AUTOPREFETCH: list variables for the file ''%s'', block=%d',o.mfilename,o.mblock)
    end
    STATE.(o.session).varlist = setdiff(o.varnames,'ans');
    if nargout, out = ~exist(o.prefetchfile,'file') || o.forcedrefresh;  end
    
elseif o.update
    % ### update() method
    % # the method calculate the diff between now and the variables listed with read()
    % # the method returns the variables names corresponding to the diff operation (note that ans is valid variable)
    validsession('update()')
    STATE.(o.session).varnames = setdiff(o.varnames,STATE.(o.session).varlist);
    if ~isempty(STATE.(o.session).varnames)
        dispf('AUTOPREFETCH: variables to be saved for file ''%s'', block=%d',o.mfilename,o.mblock)
        cprintf({'\tAUTOPREFETCH:\t[%20s]\n' '\n'},STATE.(o.session).varnames)
    else
        dispf('AUTOPREFETCH: no variables to save')
    end
    if nargout, out = STATE.(o.session).varnames; end
    
elseif o.write
    % ### write() method
    % # the method save the variables (in memory and in a prefetch file) refreshed with update()
    % # the method returns the name of the prefetch file
    validsession('write()')
    overwriting = existprefetchfile(o.prefetchfile,o.overwrite || o.forcedrefresh,'PREFETCH');
    if isempty(STATE.(o.session).varnames)
        error(['\nAUTOPREFETCH write(): nothing to write for ''%s'', block=%d.\n\tUse p.read() before calling p.write().\n' ...
               'Perhaps the script has been executed when all variables were still in memory.\n' ...
               'In this latter case, clear memory and restart the script/function.'],o.prefetchfile,o.mblock)
    elseif ~isempty(STATE.(o.session).vars) || overwriting
        dispf('\nWARNING AUTOPREFETCH: write() is overwritting previous data')
    elseif length(o.vars)~=length(STATE.(o.session).varnames)
        dispf('\nAUTOPREFETCH write() error: inconsistent number of variables\n\t--DUMP --')
        STATE.(o.session).nfo, dispf('\t--DUMP --\n')
        error('method write(), see above')
    end
    vars = cell2struct(o.vars,STATE.(o.session).varnames);
    vars.AUTOPREFETCH_NFO = STATE.(o.session).nfo;
    save(o.prefetchfile,'-struct','vars')
    rotatematfile(o.prefetchfile,o.frequency,o.nrotates,o.rotateeco)
    STATE.(o.session).vars = vars; clear vars
    dispf('AUTOPREFETCH: save prefetch file of ''%s'', block=%d (session: %s)',o.mfilename,o.mblock,o.session)
    fileinfo(o.prefetchfile)
    if nargout, out = o.prefetchfile; end
    
elseif o.save
    % ### save() method
    % # the method save all sessions as two structures STATES and MFILES
    % # the method returns the name of the saved MAT file
    o.globalfile = validprefetchfile(o.globalfile,o.path,o.overwrite,'',false);
    if isempty(STATE)
        dispf('WARNING AUTOPREFETCH: the current database is empty')
        if nargout, out = false; end
    else
        save(o.globalfile,'STATE','MFILES')
        rotatematfile(o.globalfile,o.frequency,o.nrotates,o.rotateeco)
        dispf('AUTOPREFETCH: save %d session(s)',length(fieldnames(STATE)))
        fileinfo(o.globalfile)
        if nargout, out = o.globalfile; end
    end
    
elseif o.close
    % ### close() method
    % # the method remove the curent session from memory
    % # the method returns true if the session has not been already closed
    if ~validsession('close()')
        dispf('AUTOPREFETCH WARNING: session ''%s'' has been already closed',o.session)
        if nargout, out = false; end
    else
        STATE = rmfield(STATE,o.session);
        MFILES(ismember({MFILES.session},o.session)) = [];
        if nargout, out = true; end
    end
    
    
elseif o.closeall
    % ### close() method
    % # the method close all the internal databases
    % # the method returns true if all sessions have not already been closed
    if isempty(STATE)
        disp('AUTOPREFETCH WARNING: all previous sessions have been already closed');
        if nargout, out = false; end
    else
        STATE = struct([]);
        MFILES = repmat(struct('mfilename','','mblock',[],'session',''),0,0);
        if nargout, out = true; end
    end
    
elseif o.disp
    % ### disp() method
    % # the method displays the properties/atrributes of the current session
    % # the method returns the object definition
    validsession('disp()')
    if nargout
        out = o;
    else
        cprintf('\t%s\n',{... list of properties
        sprintf('--- AUTOPREFETCH object (session=''%s'') ---',o.session)
        sprintf('%20s : %s (block=%d)','MFILE',o.mfilename,o.mblock)
        sprintf('%20s : %s','PREFETCH FILE',o.prefetchfile)
        sprintf('%20s : %s','GLOBAL FILE',o.globalfile)
        sprintf('%20s : %d','FORCED REFRESH',o.forcedrefresh)
        sprintf('%20s : %d','OVERWRITE',o.overwrite)
        sprintf('save with %d rotations every %4g day(s) (eco=%d)',o.nrotates,o.frequency,o.rotateeco)
            } )
        if ~isempty(STATE.(o.session).varnames)
            dispf('\tLISTED VARIABLES:')
            cprintf({'%-3.8s,' '\n'},STATE.(o.session).varnames)
        end
        if ~isempty(STATE.(o.session).vars)
            dispf('\t%d VARIABLES have been saved',length(fieldnames(STATE.(o.session).vars)))
        end
        if exist(o.prefetchfile,'file')
            cprintf({'\tprefetch file saved on %s' '\n'},substructarray(fileinfo(o.prefetchfile),'date'))
        end
        if exist(o.globalfile,'file')
            cprintf({'\tglobal file saved on %s' '\n'},substructarray(fileinfo(o.globalfile),'date'))
        end        
        dispf('\t--- %d sessions in the database ---',length(MFILES))
    end
elseif o.help
    % ### help() method
    % # the method lists all available properties and attributes
    listofprop = { ... table content without layout
    '''keyword''' 'type' 'description' 'default value'
    '' '' '' ''
    'mfilename' 'string' 'name of the M file' 'mfilename()'
    'mblock' 'integer' 'block index of the Mfile' 'last block +1'
    'path' 'string' 'default path' ''''''
    'prefetchfile' 'string' 'name of the prefetch file' 'AUTOPREFETCH_'
    'globalfile' 'string' 'name of the global file' 'AUTOPREFETCH_ALL'
    'frequency' 'scalar' 'frequency of file rotation in days' '1'
    'nrotates' 'integer' 'number of copies of prefetch and global files' '2'
    'rotateeco' 'boolean' 'if true only NULL files are saved at creation' '1'
    'overwrite' 'boolean' 'if true all files will be overwritten' 'false'
    'forcedrefresh' 'boolean' 'if true prefetch files are forced to be recalculated ' 'false'
    };
    colsiz = cellfun(@(s) max(cellfun(@length,s)),num2cell(listofprop,1)); % calculate the number of characters per column
    listofprop(2,:) = arrayfun(@(n) repmat('-',1,n),colsiz,'UniformOutput',false); % separator line
    fmt = sprintf('[%%%ds] [%%-%ds] [%%-%ds] [%%-%ds]\n',colsiz); % format to apply
    dispf('--- LIST OF PUBLIC PROPERTIES/ATTRIBUTES ---') % table title
    cellfun(@(s) cprintf(fmt,s),num2cell(listofprop,2)) % generate table
    dispf('--- all keywords can be used for constructors and methods ---') % table bottom

else
    % ### unidentified method
    error('AUTOPREFETCH: unrecognized method.n\tPlease use help autoprefetch or doc autoprefetch.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISVALIDSESSION returns true if the session is valid
function isvalid=validsession(where)
    ok = isfield(STATE,o.session);
    if nargout
        isvalid = ok;
    elseif ~isfield(STATE,o.session)
        error('\nAUTOPREFETCH %s error: the session ''%s'' does not exist.\n\tUse p = autoprefetch(''open'') first before calling %s'...
            ,where,o.session,where)
    end
end

%VALIDPREFETCHFILE returns a valid fullpath for prefetch or global file
function f=validprefetchfile(pfile,ppath,overwrite,mfile,isprefetch)
    if nargin<5, isprefetch = true; end
    if isprefetch, typef='PREFETCH'; else typef='GLOBAL'; end
    [p,n,e] = fileparts(pfile);
    if isempty(n)
        if isprefetch, n='commandline';
        else n=sprintf('AUTOPREFETCH_ALL_%s.mat',datestr(now,'yyyy-mm-dd_HH_MM'));
        end
    end
   if isempty(p)
       if isempty(ppath) && isprefetch, ppath =fileparts(which(mfile)); end
       if isempty(ppath), ppath =tempdir; end
       p = ppath;
   end
   if ~exist(p,'dir'), error('AUTOPREFETCH: the path to be used for %s file ''%s'' does not exist',typef,p), end
   if isempty(e), e = '.mat'; end
   f = fullfile(p,[n e]);
   existprefetchfile(f,overwrite,typef)
end

%EXISTPREFETCHFILE returns true if the file exists
% generates an error if the file exists and if overwrite is not set
function overwriting = existprefetchfile(f,overwrite,typef)
    fchk = exist(f,'file');
    if ~overwrite && fchk
        if nargin<3, typef='PREFETCH'; end
        error('AUTOPREFETCH: the following %s file already exist:\n%s\n\tSet ''overwrite'' flag to remove this message',typef,f)
    end
    if nargout, overwriting = fchk; end
end

end