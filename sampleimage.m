function [Iout,oout] = sampleimage(I,varargin)
%SAMPLEIMAGE extract and optimize an image for prez like movies
%   syntax: Isample = sampleimage(I [options,'property',value,...])
%           [Isample,options] = sampleimage(...)
%   I: very high resolution RGB image
%   options: options used for the current image
%
%   List of porperties (property:default value)
%                  x: center position along x
%                  y: center position along y
%               zoom: zoom factor (>1 to zoom)
%              width: width of the final image (use a constant value to generate movie) (default = 800)
%             height: (default = 600)
%             filter: rendering filter (default = fspecial('Gaussian',3))
%             method: resize method (default = 'bicubic')
%       Antialiasing: antialiasing flag (default = true)
%              angle: rotation angle (default = 0)
%     rotationmethod: rotation method (default = 'bicubic')
%
%   Example:
%       local = 'C:\Data\Olivier\INRA\Etudiants & visiteurs\Mai Nguyen\sandbox';
%       imfile = 'bitmap.png';
%       I = imread(fullfile(local,imfile));
%       figure, imshow(I), axis([5500 6500 5500 6500])
%       figure, imshow(sampleimage(I,'x',5900,'y',6070,'zoom',10));


%INRA\MS 2.1 - 17/05/2014 - Olivier Vitrac - rev.

% TODO LIST / KNOWN ISSUE
%   beyond the physical resolution of the image, rotation causes a random translation (not fixed at this time)



% arg check
if nargin<1, error('one image is required'), end
[H,W,d] = size(I);
%classI = class(I);

% default
default = struct(...
   'x',round(W/2),...
   'y',round(H/2),...
   'zoom',1,...
   'width',800,...
   'height',600,...
   'filter',fspecial('Gaussian',3),...
   'method','bicubic',...
   'Antialiasing',true,...
   'angle',0,...
   'rotationmethod','bicubic' ...
   );

% options
o = argcheck(varargin,default);
scale = max(o.width/W,o.height/H); % scale if the full image is plotted at with a size height x width
maxzoom = 1/scale;

% do zoom
if o.zoom<maxzoom % no zoom, just undersampling
    w = round(min(W,o.width/(o.zoom*scale)));
    h = round(min(H,o.height/(o.zoom*scale)));
    if isempty(o.angle) || o.angle ==0 % without rotation
        Icrop = imcrop(I,[o.x-round(w/2) o.y-round(h/2) w h]);
    else
        Icrop = imcrop(I,[o.x-w o.y-h 2*w 2*h]);
        Icrop = imrotate(Icrop,o.angle,o.rotationmethod,'crop'); % do rotation (if any)
        Icrop = imcrop(Icrop,[round(w/2) round(h/2) w h]);
    end
    Iout = imresize(Icrop,[o.height o.width],o.method,'Antialiasing',o.Antialiasing);
else
    realzoom = o.zoom/maxzoom;
    w = round(min(W,o.width));
    h = round(min(H,o.height));
    if isempty(o.angle) || o.angle ==0 % without rotation
        Icrop = imcrop(I,[o.x-round(w/(2*realzoom)) o.y-round(h/(2*realzoom)) w h]);
    else
        Icrop = imcrop(I,[o.x-round(w/realzoom) o.y-round(h/realzoom) 2*w 2*h]);
        Icrop = imrotate(Icrop,o.angle,o.rotationmethod,'crop'); % do rotation (if any)
        Icrop = imcrop(Icrop,[round(w/(2*realzoom)) round(h/(2*realzoom)) w h]);
    end
    Iout = imresize(Icrop,realzoom,o.method,'Antialiasing',o.Antialiasing,'OutputSize',[o.height o.width]);
end

% filtering
if ~isempty(o.filter)
    Iout = imfilter(Iout,o.filter);
end

% additional output
if nargout>1, oout = o; end