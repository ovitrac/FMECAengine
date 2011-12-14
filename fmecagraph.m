function [hgraphtmp,hparentobjout] = fmecagraph(fmecadb,values,varargin)
%FMECAGRAPH plots a graph from a FMECAdatabase (as FMECAENGINE does)
%   SYNTAX: fmecagraph(fmecadb [,values,property1,propvalue1,property2,propvalue2,...])
%           hg = fmecagraph(...) to retrieve the handle of the axes that contain the graph (not interactive)
%       fmecadb: FMECA database (output of fmeaengine)
%        values: values to set a colorscale such values.(nodes{i}) = value of the ith node
%                default = 0;
%       Implemented pair properties values
%          parent: 'parent' (default) or 'inherit'
%             min: minimal value for color scale
%             max: maximal value for color scale
%        colormap: colorscale setup between min and max (default = jet(64))
%defaulfacetcolor: face color for out of bounds values (default = [0 0 0])
%defaultedgecolor: edge color for out of bounds values (default = [1 0 0])
%defaultlinewidth: linewidth for out of bounds values (default = 2)
%  defaultmagnify: box magnification for out of bounds values (default = 1.2)
%     operation: operation to apply when several values are found (default = @(x) max(x))
%
%   OPTIONS: [hg,hbiograph] = fmecagraph(...) returns also the original biograph object (interactive)
%
%   TIP: use content = gcfd(hg); and figure, scfd(content,'noaxes','nolengend') to recopy the graph into new axes (e.g. a subplot)
%
%
%   See also: PNGTRUNCATEIM, FMECAENGINE, FMECASINGLE, GCFD, SCFD

% Migration 2.0 - 24/05/11 - INRA\Olivier Vitrac - rev. 14/12/11

% Revision history
% 14/12/11 add parent as property, update help to enable the copy of a graph

% Default
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
    'paperproperties',struct('PaperUnits','Centimeters','PaperType','A0','PaperOrientation','Landscape'));

% arg check
if nargin<1, error('1 inputs are required'), end
if nargin<2, values = struct([]); end
options = argcheck(varargin,default);
if ~isstruct(fmecadb) || numel(fmecadb)>1, error('fmecadb must be created with fmecaengine'), end
fdb = fieldnames(fmecadb); nfdb = length(fdb);
if isempty(values)
    fvalues = fdb; nvalues = numel(fvalues);
    values = cell2struct(num2cell(zeros(nvalues,1),2),fvalues,1);
else
    fvalues = fieldnames(values); nvalues = numel(fvalues);
end
missing = setdiff(fvalues,fdb);
if ~isempty(missing), cellfun(@(x) dispf('\t''%s'' is missing',x),missing),error('some fields in values are not available in fmecadb (SEE ABOVE)'); end
[~,i] = intersect(fdb,fvalues); order = false(nfdb,1); order(i) = true; fvalues = fdb(order); % order fvalues as fdb
if isinf(options.min), options.min = min(cellfun(@(x) min(values.(x)),fvalues)); end
if isinf(options.max), options.max = max(cellfun(@(x) max(values.(x)),fvalues)); end
lvalues = cellfun(@(x) length(values.(x)),fvalues);
if any(lvalues>1), dispf('WARNING: several values were found by fields'), end
vvalues = cellfun(@(x) options.operation(values.(x)),fvalues);

% color scale
col = interp1(linspace(0,1,size(options.colormap,1)),options.colormap,(vvalues-options.min)/(options.max-options.min),'cubic');
bad = (vvalues<options.min) | (vvalues>options.max);
col(bad,:) = repmat(options.defaultfacecolor,length(find(bad)),1);

% build graph
codes = cell2struct(num2cell((1:nvalues)',2),fvalues);
g = sparse(nvalues,nvalues);
for i=1:nvalues
    p = fmecadb.(fvalues{i}).(options.parent);
    if ~isempty(p), g(codes.(p),codes.(fvalues{i})) = 1; end
end

% plot
gobj = biograph(g,fvalues);
hobj=view(gobj);
for i=1:nvalues
    set(hobj.Nodes(i),'color',col(i,:))
    if sum(col(i,:))/3<.4, set(hobj.Nodes(i),'TextColor',[1 1 1]); else set(hobj.Nodes(i),'textcolor',[0 0 0]); end
    if bad(i)
        set(hobj.Nodes(i),'LineColor',options.defaultedgecolor,'Size',options.defaultmagnify*get(hobj.Nodes(i),'Size'),'LineWidth',options.defaultlinewidth)
    end
end


% output
if nargout
    hparentobj = get(hobj.hgAxes,'Parent');  set(hparentobj,options.paperproperties)
    hgraphtmp = figure('Units','Points','PaperPosition',get(hparentobj,'PaperPosition'),options.paperproperties);
    copyobj(hobj.hgAxes,hgraphtmp); set(hgraphtmp,'Units','Pixels');
end

if nargout>1, hparentobjout = hparentobj; end