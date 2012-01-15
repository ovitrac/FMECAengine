function hout = key2keygraph(keytree,hax,varargin)
%KEY2KEYGRAPH plots a key2key query as a tree (based on fmecagraph)
%   Syntax: hax = key2keygraph(keyquery [,hax,'parameter1',value1,'parameter2',value2,...])
%      keyquery: key2key object (third output of key2key)
%           hax: axes handle to plot the tree (default=[])
%
% ------------------------------------------------------------------
%   keyquery: structure with fields node1,node2...noden
% ------------------------------------------------------------------
%   nodei is a sub-structure with fields:
%              parent: parent node (string)
%                 key: subkey (string)
%              values: result of the subkey (double array or char array)
%          isterminal: flag (true if the node is terminal)
%  TIP: set isterminal flag to false to remove the display of final or intermediate query results
%
% ------------------------------------------------------------------
%   Pair property/value: (default=default_value)
%   Properties derived from fmecagraph are noted with [*]
%   see fmecagrah for further details
% ------------------------------------------------------------------
%         'alignment': HorizontalAlignement (default='center') [*]
%          'colormap': use this property to change colors      [*]
%          'fontsize': fontsize (default=10)                   [*]
%       'layoutscale': layout scale (default=1)                [*]
%        'layouttype': tree layout (default='hierarchical')    [*]
%                      alternatives: 'radial', 'equilibrium'
%           'nocolor': flag to remove color scale (default=false)
%      'placeholders': (default='')                            [*]
%             'scale': (default=1)                             [*]
%        'shapenodes': (default='box')                         [*]
%'shapeterminalnodes': (default='box')                         [*]
%         'sizenodes': (default=[])                            [*]
%            'resize': (default=[])                            [*]
% 'sizeterminalnodes': (default=[])                            [*]
%'terminalexpression': structure setting expressions to display terminal nodes
%           default.numsingle =  @(x)sprintf('1 value\n(%0.4g)',x)
%           default.nummore   =  @(x)sprintf('%d values\n(%0.4g\\leq%0.4g\\leq%0.4g)
%        default.stringsingle =  @(x)sprintf('1 value\n(%s)',x)
%        default.stringmore   =  @(x)sprintf('%d values\n(strings)',length(x))
%     'textwidth': text width for text wrap                    [*]
%       'texttol': lookaround before wrapping text             [*]
%
%   See also: KEY2KEY, FMECAGRAPH, BUILDMARKOV
%
%   Example for Mai (for her thesis)
%   figure,  ha=subplot(3,3,[2 3 5 6]); key2keygraph(detailsD,ha,'fontsize',5,'terminalexpression',struct('nummore',@(x) sprintf('%d values\n\\fontsize{6}(%0.4g\\leq%0.4g\\leq%0.4g',length(x),min(x),median(x),max(x))))

% Migration 2.0 - 24/12/2011 - INRA\Olivier Vitrac - rev. 13/01/2012

% Revision history
% 27/12/2011 release candidate
% 28/12/2011 add shapenodes, sizenodes
% 29/12/2011 improved help, fix cla
% 30/12/2011 add resize, colormap, fix figure resize on screen and paper
% 13/01/2012 fix axes, apply also argcheck to terminalexpression

% Default parameters
eol = char(10);
terminalexpression = struct(...
    'numsingle',@(x) sprintf('1 value\n(%0.4g)',x),...
    'nummore',@(x) sprintf('%d values\n\\fontsize{8}(%0.4g\\leq%0.4g\\leq%0.4g)',length(x),min(x),median(x),max(x)),...
    'stringsingle',@(x) sprintf('1 value\n(%s)',x),...
    'stringmore',@(x) sprintf('%d values\n(strings)',length(x)) ...
    );
default = struct(...
    'colormap',[],...
    'nocolor',false,...
    'fontsize',10,...
    'textwidth',20,...
    'texttol',5,...
    'alignment','center',...
    'layouttype','hierarchical',... hierarchical
    'layoutscale',1,...
    'placeholders',[],...
    'resize',1,...
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
options = argcheck(varargin,default,'','nostructexpand');
options.terminalexpression = argcheck(options.terminalexpression,terminalexpression);
nodetreelist = fieldnames(keytree)';
isterminal = cellfun(@(f) keytree.(f).isterminal, nodetreelist);

% check handle
if isempty(hax)
    hfig = gcf; hax = gca; cla, newax = true;
else
    if ~ishandle(hax) || ~strcmp(get(hax,'Type'),'axes'), error('hax must be a valid axes handle'), end
    hfig = get(hax,'Parent');
    figure(hfig), subplot(hax)
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
[hkey,hbiograph,gbiograph]= fmecagraph(keytree,numval,'names',names,...
      'placeholders',options.placeholders,...
      'terminalnodes',terminalnodes,...
      'fontsize',options.fontsize,'alignment',options.alignment,'colormap',colormap,...
      'layouttype',options.layouttype,'layoutscale',options.layoutscale,...
      'scale',options.scale,'resize',options.resize,...
      'shapenodes',options.shapenodes,'sizenodes',options.sizenodes,...
      'shapeterminalnodes',options.shapeterminalnodes,'sizeterminalnodes',options.sizeterminalnodes);
  
% Store important properties before copying data to a virgin figure/axes
set(hbiograph,'Units','pixels')
position = get(hbiograph,'Position');
delete(hbiograph) % keep only the copy
paperunits = get(hkey,'PaperUnits');
paperposition = get(hkey,'PaperPosition');

% Clear copy to hax
gkey = gcfd(hkey);
delete(hkey);
figure(hfig), subplot(hax), scfd(gkey,'noaxes','nolegend')
set(hax,'visible','off','Xlim',gkey.Xlim,'Ylim',gkey.Ylim);

% Update the position of axes and figure margins
if newax
    set(hax,'Units',gkey.Units,'Position',gkey.Position)
    set(hfig,'PaperUnits',paperunits);
    paperdefault = get(hfig,'paperposition');
    set(hfig,'Position',position,'PaperPosition',[paperdefault(1:2)-(paperposition(3:4)-paperdefault(3:4))/2,paperposition(3:4)])
end

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
