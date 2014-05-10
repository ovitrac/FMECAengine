function hout = treedir(myfolder)
%TREEDIR displays the content of a folder using uitree
%   Syntax: treedir
%           hfig = treedir('myfolder')

% MS 2.1 - 08/03/2013 - INRA\Olivier Vitrac - rev.
% Source: http://undocumentedmatlab.com/blog/uitree/

% arg check
if nargin<1, myfolder = pwd; end

hPanel = figure;
% panel0=uipanel('Parent', hPanel, 'Title', myfolder);
[~, container] = uitree('v0', 'Root',myfolder, 'Parent',hPanel); % Parent is ignored
set(container, 'Parent', hPanel,'Units','normalized','position',[0 0 1 1]);  % fix the uitree Parent
set(hPanel,'ToolBar','none')
%set(hPanel,'MenuBar','none')

% output
if nargout, hout = hPanel; end