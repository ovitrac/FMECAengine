function hout = key2keygraph(keytree,hax,varargin)
%KEY2KEYGRAPH plots a key2key query as a tree (based on fmecagraph)
%   Syntax: hax = key2keygraph(keytree [,hax,'parameter1',value1,'parameter2',value2,...])
%      keytree: key2key object (see third output of key2key)
%           hax: axes handle to plot the tree (default=[])
%          
%   Pair property/value: (default=default_value)
%   (properties derived from fmecagraph are noted with [*])
%         'alignment': HorizontalAlignement (default='center') [*]
%          'fontsize': fontsize (default=10)                   [*]
%       'layoutscale': layout scale (default=1)                [*]
%        'layouttype': tree layout (default='hierarchical')    [*]
%           'nocolor': flag to remove color scale (default=false)
%      'placeholders': (default='')                            [*]
%             'scale': (default=1)                             [*]
%        'shapenodes': (default='box')                         [*]
%'shapeterminalnodes': (default='box')                         [*]
%         'sizenodes': (default=[])                            [*]
% 'sizeterminalnodes': (default=[])                            [*]
%'terminalexpression': structure setting expressions to display terminal nodes
%           default.numsingle =  @(x)sprintf('1 value\n(%0.4g)',x)
%           default.nummore   =  @(x)sprintf('%d values\n(%0.4g\\leq%0.4g\\leq%0.4g)
%        default.stringsingle =  @(x)sprintf('1 value\n(%s)',x)
%        default.stringmore   =  @(x)sprintf('%d values\n(strings)',length(x))
%     'textwidth': text width for text wrap                    [*]
%       'texttol': lookaround before wrapping text             [*]

%   See also: KEY2KEY, FMECAGRAPH

% Migration 2.0 - 06/12/2011 - INRA\Olivier Vitrac - rev. 

% Default parameters
eol = char(10);
terminalexpression = struct(...
    'numsingle',@(x) sprintf('1 value\n(%0.4g)',x),...
    'nummore',@(x) sprintf('%d values\n\\fontsize{8}(%0.4g\\leq%0.4g\\leq%0.4g)',length(x),min(x),median(x),max(x)),...
    'stringsingle',@(x) sprintf('1 value\n(%s)',x),...
    'stringmore',@(x) sprintf('%d values\n(strings)',length(x)) ...
    );
default = struct(...
    'nocolor',false,...
    'fontsize',10,...
    'textwidth',20,...
    'texttol',5,...
    'alignment','center',...
    'layouttype','hierarchical',... hierarchical
    'layoutscale',1,...
    'placeholders',[],...
    'scale',1,...
    'shapenodes','box',...
    'shapeterminalnodes','invtrapezium',...
    'sizenodes',[],...
    'sizeterminalnodes',[],...
    'terminalexpression',terminalexpression ...
    );

% arg check
if nargin<1, error('one argument is at least required'), end
if nargin<2, hax = []; end
if ~isstruct(keytree), error('keytree must be a structure'); end
if ~all(cellfun(@(f) all(isfield(keytree.(f),{'parent' 'key'})),fieldnames(keytree)))
    error('keytree must be created by key2key (3rd output), see KEY2KEY ')
end
options = argcheck(varargin,default);
nodetreelist = fieldnames(keytree)';
isterminal = cellfun(@(f) keytree.(f).isterminal, nodetreelist);

% check handle
if isempty(hax)
    hfig = gcf; hax = gca; newax = true;
else
    if ~ishandle(hax) || strcmp(get(hax,'Type'),'axes'), error('hax must be a valid axes handle'), end
    hfig = get(hax,'Parent');
    figure(hfig), subplot(hax)
    cla
    newax = false;
end

% wrap text
names = cell2struct(...
    cellfun(@(f) wraptext(keytree.(f).key,options.textwidth,eol,options.texttol,true),nodetreelist,'UniformOutput',false),...
    nodetreelist,2 );

% values
numval = cell2struct(cellfun(@(f) length(keytree.(f).values),nodetreelist,'UniformOutput',false),nodetreelist,2);

% terminal nodes
terminalnodes = rmfield(numval,nodetreelist(~isterminal));
for f=nodetreelist(isterminal)
    val = keytree.(f{1}).values;
    if ischar(val), val = {val}; end
    if iscellstr(val)
        if length(val)<2
            terminalnodes.(f{1}) = options.terminalexpression.stringsingle(val{1});
        else
            terminalnodes.(f{1}) = options.terminalexpression.stringmore(val);
        end        
    else
        if length(val)<2
            terminalnodes.(f{1}) = options.terminalexpression.numsingle(val);
        else
            terminalnodes.(f{1}) = options.terminalexpression.nummore(val);
        end
    end
end

% plot tree
if options.nocolor, numval = []; end
[hkey,hbiograph]= fmecagraph(keytree,numval,'names',names,'placeholders',options.placeholders,'terminalnodes',terminalnodes,...
      'fontsize',options.fontsize,'alignment',options.alignment,...
      'layouttype',options.layouttype,'layoutscale',options.layoutscale,'scale',options.scale,...
      'shapenodes',options.shapenodes,'sizenodes',options.sizenodes,'shapeterminalnodes',options.shapeterminalnodes,'sizeterminalnodes',options.sizeterminalnodes);
delete(hbiograph) % keep only the copy

% Copy again the copy to hax
gkey = gcfd(hkey); delete(hkey);
figure(hfig), subplot(hax), scfd(gkey,'noaxes','nolegend')
set(hax,'visible','off','Xlim',gkey.Xlim,'Ylim',gkey.Ylim);
if newax, set(hax,'Units',gkey.Units,'Position',gkey.Position), end

% output
if nargout, hout = hax; end

%%%%%%%%%%%% OBSOLETE
% % layoutoverride
% if isempty(options.sizenodes)
%     nrows = max(cellfun(@(f) length(find(names.(f)==eol))+1,nodetreelist));
%     widths = max(cellfun(@(f) max(cellfun('length',regexp(names.(f),'\n','split'))),nodetreelist));
%     if ~isempty(options.fontsize)
%         options.sizenodes = [widths*12 nrows*2];
%     else
%         options.sizenodes = [widths nrows]*options.fontsize;
%     end
% end

% placeholders
% placeholders = cell2struct(repmat({repmat('#',1,options.textwidth+options.texttol)},1,nnodes),nodetreelist,2);
