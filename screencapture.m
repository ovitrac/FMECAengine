function out=screencapture(varargin)
%SCREENCAPTURE capture the content of a region, a figure or of the screen
%       syntax: ut=screencapture('parameter1',value1,'parameter1',value1,keyword,...)
%       Pair property/value
%           position: [x y width height] (default = screen size)
%             figure: figure handle (default='')
%           filename: full filename (default = /tmp/screencapture.png)
%                bar: height of the bar title (default = 22)
%      cameratoolbar: height of the cameratoolbar (default = 27)
%       Keyword
%           'truncate' to apply pngtruncateim to the capture
%
%   Note: there is a pure Matlab alternative but usually less efficient
%       tmp = getframe; imwrite(tmp.cdata,'myfile','png')


% Migration 2.1 - 05/02/12 - INRA\Olivier Vitrac - rev.

% default
default = struct(...
    'position',get(0,'ScreenSize'),...
    'figure',[],...
    'filename',fullfile(tempdir,'screencapture.png'),...
    'bar',22,...
    'cameratoolbar',27 ...
    );
keywords = 'truncate';

% argcheck
options = argcheck(varargin,default,keywords);
if ~isempty(options.figure)
    options.figure = options.figure(1);
    if ~ishandle(options.figure) || ~strcmp(get(options.figure,'type'),'figure'), error('a valid figure handle is required for figure'); end
    menu = get(options.figure,'menu');
    toolbar = get(options.figure,'toolbar');
    cameratoolbarvis = cameratoolbar('GetVisible');
    set(options.figure,'menu','none','toolbar','none')
    if cameratoolbarvis, cameratoolbar('Hide'), options.bar = options.bar+options.cameratoolbar; end
    set(options.figure,'Units','pixels')
    figure(options.figure), drawnow
    options.position = get(options.figure,'position');
    options.position(2) = options.position(2) + options.bar;
    options.position(4) = options.position(4) - options.bar;
end
[~,~,ext] = fileparts(options.filename);
if isempty(ext), options.filename = [options.filename '.png']; end

% snapshot derived from:
% http://www.mathworks.com/matlabcentral/fileexchange/11363-screencapture
% http://www.public.iastate.edu/~java/docs/api/java.awt.Rectangle.html
area = java.awt.Rectangle(options.position(1)-1,options.position(2)-1,options.position(3),options.position(4));
robot = java.awt.Robot;
image = robot.createScreenCapture(area);
filehandle = java.io.File(options.filename);
javax.imageio.ImageIO.write(image,'png',filehandle)
if options.truncate, pngtruncateim(options.filename,0,0); end
out = options.filename;

% restore figure properties
if ~isempty(options.figure)
    set(options.figure,'menu',menu,'toolbar',toolbar)
     if cameratoolbarvis, cameratoolbar('Show'), end
end

%Initial code
% robo = java.awt.Robot;
% t = java.awt.Toolkit.getDefaultToolkit();
% rectangle = java.awt.Rectangle(t.getScreenSize());
% image = robo.createScreenCapture(rectangle);
% filehandle = java.io.File('screencapture.jpg');
% javax.imageio.ImageIO.write(image,'jpg',filehandle);
% imageview('screencapture.jpg');
%
% rasta=image.getRGB(0,0,w,h,[],0,w); %get RGB data from bufferedimage 
% %convert java color integers to matlab RGB format:
% rasta=256^3+rasta;
% B=uint8(mod(rasta,256));
% G=uint8(mod((rasta-int32(B))./256,256));
% R=uint8(mod((rasta-256*int32(G))./65536,256));