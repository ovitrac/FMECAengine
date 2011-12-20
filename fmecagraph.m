function [hgraphtmp,hparentobjout] = fmecagraph(fmecadb,values,varargin)
%FMECAGRAPH plots a graph from a FMECAdatabase (as FMECAENGINE does)
%   SYNTAX: fmecagraph(fmecadb [,values,property1,propvalue1,property2,propvalue2,...])
%           hg = fmecagraph(...) to retrieve the handle of the axes that contain the graph (not interactive)
%         fmecadb: FMECA database (output of fmeaengine)
%          values: values to set a colorscale such values.(nodes{i}) = value of the ith node
%                default = 0;
%       Implemented pair properties values
%            parent: 'parent' (default) or 'inherit'
%               min: minimal value for color scale
%               max: maximal value for color scale
%          colormap: colorscale setup between min and max (default = jet(64))
%  defaulfacetcolor: face color for out of bounds values (default = [0 0 0])
%  defaultedgecolor: edge color for out of bounds values (default = [1 0 0])
%  defaultlinewidth: linewidth for out of bounds values (default = 2)
%    defaultmagnify: box magnification for out of bounds values (default = 1.2)
%         operation: operation to apply when several values are found (default = @(x) max(x))
%             names: alternate node names for plotting coded as names.(nodes{i})='alternate name'
%                    used either in orginal (interactive) graph and copied graph
%      placeholders: as names but to force a precribed design
%     terminalnodes: terminal values plotted as terminal nodes codes as terminalnodes.(nodes{i})='some text'
%           weights: value connecting between nodes{i} and nodes{i+1} (default=1) coded as weights.(nodes{i})=value
%                    'auto' forces weights to be calculated as: values.(nodes{i})-values.(nodes{i-1})
%shapeterminalnodes: shape for terminal nodes (default = 'ellipse')
%         rootvalue: root value to be used for calculating weights
%          fontsize: as text()
%         
%
%   OPTIONS: [hg,hbiograph] = fmecagraph(...) returns also the original biograph object (interactive)
%
%   TIP: use content = gcfd(hg); and figure, scfd(content,'noaxes','nolengend') to recopy the graph into new axes (e.g. a subplot)
%
%
%   See also: PNGTRUNCATEIM, FMECAENGINE, FMECASINGLE, GCFD, SCFD

% Migration 2.0 - 24/05/11 - INRA\Olivier Vitrac - rev. 19/12/11

% Revision history
% 14/12/11 add parent as property, update help to enable the copy of a graph
% 18/12/11 fix NaN as color, add names, weights, terminalnodes, shapeterminalnodes
% 19/12/11 fix weights calculations, add placeholders, fontsize

% Default
autoweights = false;
default = struct(...
    'parent','parent',...
    'min',-Inf,...
    'max',Inf,...
    'colormap',jet(64),...
    'defaultfacecolor',[0 0 0],...
    'defaultedgecolor',[1 0 0],...
    'defaultlinewidth',2,...
    'defaultmagnify',1.2,...
    'operation',@(x) max(x),...
    'paperproperties',struct('PaperUnits','Centimeters','PaperType','A0','PaperOrientation','Landscape'),...
    'names',struct([]),...
    'placeholders',struct([]),...
    'weights',[],...
    'terminalnodes',[],...
    'shapeterminalnodes','ellipse',...
    'rootvalue',0,...
    'fontsize',10);
minplaceholderlength = 4;
placeholderidx = sprintf('%%0.%dd',minplaceholderlength-1);

% arg check
if nargin<1, error('1 inputs are required'), end
if nargin<2, values = struct([]); end
options = argcheck(varargin,default,'','nostructexpand');
if ~isstruct(fmecadb) || numel(fmecadb)>1, error('fmecadb must be created with fmecaengine'), end
fdb = fieldnames(fmecadb); nfdb = length(fdb);
if isempty(options.names) && ~isstruct(options.names), options.names = struct([]); end
if isempty(options.placeholders) && ~isstruct(options.placeholders), options.placeholders = struct([]); end
if isempty(options.terminalnodes) && ~isstruct(options.terminalnodes), options.terminalnodes= struct([]); end
if ischar(options.weights) && strcmpi(options.weights,'auto'), autoweights = true; options.weights = []; end
if isempty(options.weights), options.weights = cell2struct(num2cell(ones(nfdb,1),2),fdb); end
if ~isstruct(options.names), error('names must be a structure such as names.step = ''alternate name'''), end
if ~isstruct(options.placeholders), error('placeholders must be a structure such as names.step = ''some text to fill space'''), end
if ~isstruct(options.terminalnodes), error('names must be a structure such as terminalnodes.step = ''some text'''), end
if ~isstruct(options.weights), error('names must be a structure such as weights.step = value'), end

% extract values to match fdb
if isempty(values)
    fvalues = fdb; nvalues = numel(fvalues);
    values = cell2struct(num2cell(zeros(nvalues,1),2),fvalues,1);
else
    for f=fdb', if ~isfield(values,f{1}), values.(f{1}) = NaN; end, end
    values = orderfields(values,fdb); fvalues = fieldnames(values); nvalues = numel(fvalues);
    if ~cellcmp(fdb,fvalues), error('fmecadb and values must share the same fields'); end
end
if isinf(options.min), options.min = min(cellfun(@(x) min(values.(x)),fvalues)); end
if isinf(options.max), options.max = max(cellfun(@(x) max(values.(x)),fvalues)); end
lvalues = cellfun(@(x) length(values.(x)),fvalues);
if any(lvalues>1), dispf('WARNING: several values were found by fields'), end
vvalues = cellfun(@(x) options.operation(values.(x)),fvalues);

% color scale
col = interp1(linspace(0,1,size(options.colormap,1)),options.colormap,(vvalues-options.min)/(options.max-options.min),'cubic');
bad = (vvalues<options.min) | (vvalues>options.max);
col(bad,:) = repmat(options.defaultfacecolor,length(find(bad)),1);

% extract weights to match fdb
for f=fdb', if ~isfield(options.weights,f{1}), options.weights.(f{1}) = 1; end, end
options.weights = orderfields(options.weights,fdb);
if ~cellcmp(fdb,fieldnames(options.weights)), error('fmecadb and options.weights must share the same fields'); end
w = cellfun(@(f) options.weights.(f)(1),fdb);

% terminal nodes
fdbterminal = fdb(cellfun(@(f) fmecadb.(f).isterminal,fdb));
unknown = setdiff(fieldnames(options.terminalnodes),fdbterminal);
if ~isempty(unknown)
    dispf('ERROR\t%d prescribed terminalnodes are not terminal',length(unknown))
    cellfun(@(f) dispf('\t''%s'' is not terminal',f),unknown)
    error('bad terminal or invalid nodes in terminalnodes, see above.')
end
terminalnodes = fieldnames(options.terminalnodes);
iterminalnodes = cellfun(@(f) find(ismember(fdb,f)),terminalnodes);
nterminalnodestoadd = length(terminalnodes);

% graph topology
parents = cellfun(@(f) fmecadb.(f).(options.parent),fdb,'UniformOutput',false);
[~,~,m,c] = buildmarkov(fdb,parents); % m=map, c= list of children

% auto weights: FMECAengine method (path-by-path)
if autoweights
    for jpath=1:size(m,2)
        isteps = m(m(:,jpath)>0,jpath);
        w(isteps)=diff([options.rootvalue;vvalues(isteps)]);
    end
    w(w==0) = NaN;
end

% build graph
g      = sparse(nvalues+nterminalnodestoadd,nvalues+nterminalnodestoadd); %g = sparse(nvalues,nvalues);
[gnames,glabels,glabelscopy] = deal(cell(nvalues+nterminalnodestoadd,1));
gnames(1:nvalues) = fvalues;
glabels(1:nvalues) = fvalues;
glabelscopy(1:nvalues) = fvalues;
% obsolete method (value-by-value)
% codes = cell2struct(num2cell((1:nvalues)',2),fvalues);
% for i=1:nvalues
%     p = fmecadb.(fvalues{i}).(options.parent);
%     if ~isempty(p), g(codes.(p),codes.(fvalues{i})) = 1; end
% end
% 
% main nodes: as FMECAengine does
for i=1:nvalues
    j=c(:,i)>0; g(i,c(j,i))=w(i); %#ok<SPRIX>
    if isfield(options.placeholders,fvalues{i})
        glabels{i} = options.placeholders.(fvalues{i});
        if length(glabels{i}) < minplaceholderlength, glabels{i} = [glabels{i} repmat('#',1,minplaceholderlength-length(glabels{i}))]; end
        glabels{i}(1:minplaceholderlength-1) = sprintf(placeholderidx,i);
        if isfield(options.names,fvalues{i}) % main node with alternate name
            glabelscopy{i} = options.names.(fvalues{i});
        end
    else
        if isfield(options.names,fvalues{i}) % main node with alternate name
            glabels{i} = options.names.(fvalues{i});
            glabelscopy{i} = glabels{i};
        end
    end
end

% terminal nodes
if nterminalnodestoadd
    virtualnodesnames = arrayfun(@(i) sprintf('VirtualNode%d',i),(1:nterminalnodestoadd)','UniformOutput',false);
    for i=1:nterminalnodestoadd
        g(iterminalnodes(i),nvalues+i) = w(iterminalnodes(i)); %#ok<SPRIX>
        gnames{nvalues+i} = virtualnodesnames{i};
        glabels{nvalues+i} = options.terminalnodes.(terminalnodes{i});
        glabelscopy{nvalues+i} = glabels{nvalues+i};
    end
end

% create graph and apply labels
names = cell2struct(glabels,gnames);
gobj = biograph(g,gnames);
for i=1:length(gobj.Nodes)
    gobj.Nodes(i).Label = names.(gobj.Nodes(i).ID);
    gobj.Nodes(i).UserData = gobj.Nodes(i).ID;
end


% plot graph and apply layout
hobj=view(gobj);
for i=1:nvalues
    if ~isnan(col(i,:))
        set(hobj.Nodes(i),'color',col(i,:));
        if sum(col(i,:))/3<.4, set(hobj.Nodes(i),'TextColor',[1 1 1]);
        else set(hobj.Nodes(i),'textcolor',[0 0 0]);
        end
        if bad(i)
            set(hobj.Nodes(i),'LineColor',options.defaultedgecolor,'Size',options.defaultmagnify*get(hobj.Nodes(i),'Size'),'LineWidth',options.defaultlinewidth)
        end
    end
end
if nterminalnodestoadd>0, set(hobj.Nodes((1:nterminalnodestoadd)+nvalues),'Shape',options.shapeterminalnodes); end

% display edges/weights
if ~all(w==1)
    set(hobj,'ShowWeights','on')
    for i=1:length(hobj.Edges) % replace NaN by 0
        if isnan(hobj.Edges(i).Weight), hobj.Edges(i).Weight = 0; end
    end
end


% output
if nargout
    hparentobj = get(hobj.hgAxes,'Parent');  set(hparentobj,options.paperproperties)
    hgraphtmp = figure('Units','Points','PaperPosition',get(hparentobj,'PaperPosition'),options.paperproperties);
    copyobj(hobj.hgAxes,hgraphtmp); set(hgraphtmp,'Units','Pixels');    
    % replace placeholders in graph copy
    haxcopy = get(hgraphtmp,'Children');
    if ~strcmp(get(haxcopy,'Type'),'axes'), error('unable to identify axes in copied figure'), end
    hchildren = get(haxcopy,'Children');
    delete(hchildren(strcmp(get(hchildren,'visible'),'off'))) % remove hidden objects
    hchildren = get(haxcopy,'Children');
    istext = strcmp(get(hchildren,'Type'),'text');
    set(hchildren(istext),'interpreter','Tex','FontSize',options.fontsize)
    istext(istext) = cellfun(@(s) size(s,1)==1,get(hchildren(istext),'String'));
    istext(istext) = ~cellfun('isempty',regexp(strtrim(get(hchildren(istext),'String')),'\d+#+$'));
    for i=find(istext)'
        set(hchildren(i),'String',glabelscopy(str2double(regexp(get(hchildren(i),'String'),'\d+','match'))))
    end   
end

if nargout>1, hparentobjout = hparentobj; end