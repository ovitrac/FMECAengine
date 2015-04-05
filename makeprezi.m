function [s,out] = makeprezi(I,varargin)
% MAKEPREZI makes an animation from a high resolution image using continuous zooms and rotations
%   Syntax:  s = makeprezi(I [,property,value,keyword]);
%            [s,out] = makeprezi(....);
%        I: very high resolution RGB image
%      out: structure giving positions at prescribed times
%       pair/property value
%                t: 1xnt times with values between 0 (begining) and 1 (final) (default=linspace(0,1,20))
%               xy: 2xnp trajectory or anonymous function giving a 2x1 a vector
%             zoom: 1xnk kinetics or anonymous function giving a scalar (default =  @(t) 2)
%            angle: 1xna kinetics or anonymous function giving a scalar (default = @(t) 0)
%            width: image width (default=1024)
%           height: image height (default=768)
%     outputfolder: output folder (default = 'MAKEPREZI_dd-mmm-yy_HH-MM')
%      filepattern: file pattern (default = 'frame%03d')
%       Use the keyword 'print' to generate all images as files
%       Use the keyword 'inputs' to generate inputs
%
%   Example
%         local = 'C:\Data\Olivier\INRA\Etudiants & visiteurs\Mai Nguyen\sandbox';
%         imfile = 'bitmap.png';
%         I = imread(fullfile(local,imfile));
%         xy = [
%         0.2430    0.2263    0.1518    0.1135    0.1291    0.1747    0.2407    0.3367    0.4266    0.5551    0.7519    1.0098    1.1515    1.2271    1.2187    1.1923    1.0255    0.7302    0.4999    0.4386
%         0.5802    0.6270    0.6222    0.5754    0.5262    0.5106    0.5346    0.6174    0.6978    0.7458    0.7806    0.7866    0.7194    0.5874    0.4902    0.4110    0.3558    0.3318    0.2586    0.1926
%         ]*1e4; 
%         z = [0  .1 .3 .5  .9  1
%             10 10  1  20 40  20]
%         s = makeprezi(I,'t',linspace(0,1,25),'xy',xy,'zoom',z)  
%         movie(s,1,5) % movie(frames,NumberOfRepetitions,FrameRate)

% INRA\MS 2.1 - 21/05/2014 - Olivier Vitrac - rev. 03/06/2014

% Revision history
% 23/05/2014 returns out
% 23/05/2014 add height and width
% 24/05/2014 fix round xy
% 02/06/2014 prefetch s object to prevent memory fragmentation
% 03/06/2014 add filename to out

% Default
default = struct(...
    't',linspace(0,1,20)',...
    'xy', [],...
    'zoom',@(t) 2,...
    'angle',@(t) 1,...
    'outputfolder',sprintf('MAKEPREZI_%s',datestr(now,'dd-mmm-yy_HH-MM')),...
    'filepattern','frame%03d',...
    'height',768,...
    'width',1024 ...
);
keyword = {'print' 'inputs'};

% arg check
o = argcheck(varargin,default,keyword);
if ~isnumeric(o.t), error('t must be numeric'), end
if min(o.t)<0 || max(o.t)>1, error('''t'' values must be comprised between 0 and 1'), end
if isnumeric(o.xy) && size(o.xy,1)==2
    dispf('MAKEPREZI: generates a smooth trajectory from 2x%d points', size(o.xy,2))
    df = diff(o.xy,1,2);
    t = cumsum([0, sqrt([1 1]*(df.*df))]);
    t = t/t(end);
    cv = csapi(t,o.xy);
    o.xy0 = o.xy;
    o.xy = @(t) fnval(cv,t);
    o.txy = t;
    o.cvxy = cv;
else
    o.xy0 = [];
    o.txy = [];
    o.cvxy = [];
end
if ~isempty(o.zoom) && isnumeric(o.zoom)
    if length(o.zoom)==1 && ~isempty(o.txy), o.zoom = [o.txy;o.txy*0+o.zoom]; end
    if size(o.zoom,1)==2
        dispf('MAKEPREZI: generates a zoom kinetics from %d points', size(o.zoom,2))
        z = o.zoom;
        o.zoom = @(t) interp1(z(1,:),z(2,:),t,'pchip');
        o.tzoom = z(1,:);
        o.zoom0 = z(2,:);
    else
        o.tzoom = [];
        o.zoom0 = [];
    end
end
if ~isempty(o.angle) && isnumeric(o.angle)
    if length(o.angle)==1 && ~isempty(o.txy), o.angle = [o.txy;o.txy*0+o.angle]; end
    if size(o.angle,1)==2
        dispf('MAKEPREZI: generates a angle kinetics from %d points', size(o.angle,2))
        a = o.angle;
        o.angle = @(t) interp1(a(1,:),a(2,:),t,'pchip');
        o.tangle = a(1,:);
        o.angle0 = a(2,:);
    else
        o.tangle = [];
        o.angle0 = [];
    end
end
if ~isa(o.xy,'function_handle'), error('the value of ''xy'' must be a function handle'), end
if ~isa(o.zoom,'function_handle'), error('the value of ''zoom'' must be a function handle'), end
if ~isa(o.angle,'function_handle'), error('the value of ''angle'' must be a function handle'), end
if numel(o.xy(0))~=2 && numel(o.xy(1))~=2, error('''xy'' must be an anomyous value returning a 2x1 row vector for a scalar t value'), end
if numel(o.zoom(0))~=1 && numel(o.zoom(1))~=1, error('''zoom'' must be an anomyous value returning a scalar for a scalar t value'), end
if numel(o.angle(0))~=1 && numel(o.angle(1))~=1, error('''angle'' must be an anomyous value returning a scalar for a scalar t value'), end

% do animation
if ~o.inputs
    if o.print
        if exist(o.outputfolder,'dir'), error('the outputfolder already exist: ''%s''', o.outputfolder), end
        mkdir(o.outputfolder)
        dispf('MAKEPREZI: the current outpufolder is\nt\t%s',o.outputfolder)
%         fileinfo(o.outputfolder)
    end
    nt = length(o.t);
    s = repmat(struct('cdata',zeros(o.height,o.width,3,'uint8'),'colormap',[]),nt,1);
    pos = repmat(struct('t',[],'xy',[],'zoom',[],'angle',[],'filename',''),nt,1);
    screen = '';
    for i=1:nt
        t = o.t(i);     % time
        xy = round(o.xy(t)); % position
        z = o.zoom(t);  % zoom
        a = o.angle(t); % angle
        screen = dispb(screen,'MAKEPREZI: frame %d/%d with zoom=%0.4g and rotation=%0.3g',i,nt,z,a);
        s(i).cdata = sampleimage(I,'x',xy(1),'y',xy(2),'zoom',z,'angle',a,'height',o.height,'width',o.width);
        pos(i) = struct('t',t,'xy',xy,'zoom',z,'angle',a,'filename','');
        if o.print
            filename = fullfile(o.outputfolder,[sprintf(o.filepattern,i) '.png']);
            imwrite(s(i).cdata,filename)
            pos(i).filename = filename;
        end
    end
    % output
    if nargout>1, out = o; out.positions = pos; end
else
   s = o;
end