function filestr = obj2im(ho,res,filename,printoptions)
%  OBJ2IM converts any matlab object into a bitmap image which preserves the initial figure size and resolution
%    Syntax: f = obj2im([ho,resolution,filename,printoptions])
%    Inputs
%       ho: a valid handle object (default = gco), a vector of handles with a same parent, keyword
%           implemented objects include: figure, axes, graph2D, graph3D, formatted text...
%         >>  for figures, only the first one is converted
%         >>  several axes are converted in separated files
%         >>  ho='axes' is a shortcut for all axes in the current figure
%         >>  ho = 'root'generate an error
%       res: resolution in dpi (default = 300) or one of the accepted keywords
%            accepted types for res: scalar, vector, cell array
%         >> res>0 saves as image without cropping
%         >> res<0 removes the white space border (the image may be larger), threshold is applied
%         >> res='fig' saves the object as a figure instead of a bitmap
%         >> res='figxxx', where xxx = resolution value
%            saves only the content of axes as a bitmap (usefull for colormaps)
%            the axes properties are stored in data.scale
%            Note that the labels (xlabel... and title) are lost
%       filename: name of the image (default = [obj2im_ datestr(now)])
%         >> cell array are accepted for multiple axes
%       printoptions: other print options (e.g. '-opengl'), see PRINT or PRINT_JPG for details
%    Output
%       f: returned full filename (the suffix '.001', '.002'... is added when the file already exists)
%         >> f is a cell array for multiple axes
%         >> bitmap are stored in a MAT file version 7 (not compatible with prior versions)
%          The bitmap is retrieved with: data = load(f)
%               data.im: mxnxdepth RVB data of the bitmap (use image(data.im) to display it)
%                        image file
%               data.sizeim = [m n depth] size of the image
%               data.generator = 'obj2im'
%               data.resolution = resolution or 'fig'
%               data.position = position (in normalized units) of the first object copied
%               data.outerposition = outerposition (when exists, NaN if not) of the first object copied
%               data.scale = structure which store the main properties of axes, implemented properties include;
%               'xtick','ytick','ztick','xlim','ylim','zlim','xminortick','yminortick','zminortick','xticklabel'
%               'yticklabel','zticklabel','alim','xdir','ydir','zdir','xaxislocation','yaxislocation','zaxislocation'
%               'box', 'xscale','yscale'
%               data.date = datestr(now)
%
%    Example: save all axes of the current figure as 200 dpi bitmaps and recreate a similar figure consisting in several bitmaps
%       files = obj2im('axes',200,'obj2im_test') % all files are named obj2im_test.mat, obj2im_test.001.mat, obj2im_test.002.mat, ...
%       plotobj(files) % a new figure is created
%
%   See also: PLOTOBJ (to display in any axes the corresponding bitmaps)

% MS 2.0 - 13/08/07 - Olivier Vitrac - rev. 12/09/07

% Revision history
% 14/08/07 add objects of figure type (instead of an error)
% 28/08/07 add copy of several axes, update filenames instead of overwrtting them, add shortcut 'axes'
% 03/09/07 add res='fig'
% 04/09/07 force the exact copy of all labels positions, add res='figxxxx'
% 11/09/09 fix res=scalar, very large images
% 12/09/09 add sizeim

% defaults
res_default = 300;
ext = '.mat';
figext = '.fig';
tmp    = 'tmp_';
tmpfig = 'fig_';
tmpext = '.jpg';
prefix = 'obj2im_';
suffixoverwritting = 1;
suffixoverwrittingmax = 1000;
depthmax = 100;
margin = .2; % safety margin
labels = {'xlabel','ylabel','zlabel','title'}; % labels (list of recognized labels)
axesscaleproperties = {'xtick','ytick','ztick','xlim','ylim','zlim','xminortick',...
    'yminortick','zminortick','xticklabel','yticklabel','zticklabel','alim',...
    'xdir','ydir','zdir','xaxislocation','yaxislocation','zaxislocation','box',...
    'xscale','yscale'};

% arg check
if nargin<1, ho = gco; end
if ischar(ho) && strcmp(ho,'axes'), ho = flipud(findobj(get(gcf,'children'),'flat','type','axes')); end % shortcut
if nargin<2, res = []; end
if nargin<3, filename = ''; end
if nargin<4, printoptions = ''; end
if isempty(res), res=res_default; end
if ischar(res), res = {res}; end
if ~iscell(res), res = {res}; end

% check handles
if ~all(ishandle(ho)), error('invalid handle(s)'), end
o = handle(ho); no = length(o);
if no<1, error('no object to convert'), end
commonparent = double(o(1).parent);
if strcmp(o(1).type,'axes'), naxes = 1; else naxes = 0; end;
for i=2:no
    if strcmp(o(i).type,'root'), error('root cannot be copied'), end
    if strcmp(o(i).type,'axes'), naxes=naxes+1; end
    if double(o(i).parent)~=commonparent, error('all objects must have the same parent, provide the axes instead'), end
end
if naxes>1 % copy several axes
    if length(o)>naxes, error('all handles should be axes or remove the axes'), end
    filestr = cell(1,naxes);
    if ~iscell(filename), filename = {filename}; end
    filename = filename(mod((1:naxes)-1,length(filename))+1); % circular permutation if some filenames are missing
    res = res(mod(0:naxes-1,length(res))+1);
    for i=1:naxes
       filestr{i} = obj2im(ho(i),res(i),filename{i},printoptions); % partial recursion for efficiency
    end
    return
end

% check  saveasbitmap
resolution = res{1};
saveasbitmap = 1; % the entire axes is saved as bitmap
if iscell(res)
    if ischar(res{1})
        if strcmp(res{1},'fig')
            saveasbitmap = 0;
        elseif length(res{1})>3 && strcmp(res{1}(1:3),'fig')
            saveasbitmap = 2; % only the axes content is saved as bitmap
            res{1} = str2double(res{1}(4:end)); 
        else
            error('invalid resolution parameter')
        end
    end
    res = res{1};
end

% check filenames
if isempty(filename)
    cpu=cputime; cputxt=num2str((cpu-floor(cpu)));
    filename = [prefix datestr(now) '_' cputxt(3:end)]; filename(filename==':')='-';
end
[pathstr,name]= fileparts(filename);
if isempty(pathstr), pathstr = cd; end
filename = fullfile(pathstr,[name ext]);
if exist(filename,'file')
    overwritting = true;
    filenameold = filename;
    [p,f,e] = fileparts(filename);
    while suffixoverwritting<suffixoverwrittingmax && overwritting
        suffixoverwrittingstr = sprintf('.%s%d',repmat('0',1,ceil(log10(suffixoverwrittingmax))-ceil(log10((suffixoverwritting)+1))),...
            suffixoverwritting);
        filename = fullfile(p,[f suffixoverwrittingstr e]);
        suffixoverwritting = suffixoverwritting+1;
        overwritting = exist(filename,'file');
    end
    if overwritting
        error('the filename ''%s'' has been used more than %d times, clean the files by hands',filenameold,suffixoverwrittingmax)
    else
        disp(sprintf('OBJ2IM: to prevent the file ''%s'' from being overwritten,\n the suffix ''%s'' has been added',filenameold,suffixoverwrittingstr))
    end
end

% temp filename
[pathstr,name] = fileparts(filename);
if saveasbitmap
    filenametmp = fullfile(pathstr,[tmp name tmpext]);
else
    filenametmp = fullfile(pathstr,[tmpfig name figext]);
end

% check white spaces
if saveasbitmap
    removewhitespace = res<0;
    res = abs(res);
end

% object properties, copy, paste, save
if strcmp(o(1).type,'figure') % the simplest case, it is a figure, it is printed
    figure(double(o(1))) % only the first figure is retrieved
    if saveasbitmap
        print_jpg(res,filenametmp,printoptions)
    else
        saveas(double(o(1)),filenametmp,'fig')
    end
    newfig = [];
    position = o(1).position;
    outerposition = o(1).outerposition;
    scale = getprop(double(o(1)),axesscaleproperties);
else % for other cases, the objects must be extracted and copied into a new figure before printing
    % look for theirs parent axes and figure
    notfound = true;
    parent = o(1);
    ax = []; fig = [];
    depth = 1;
    while notfound && depth<depthmax
        if strcmp(parent.type,'axes'), ax = parent; end
        if strcmp(parent.type,'figure'), fig = parent; notfound = false; end
        parent = handle(parent.parent);
        depth = depth+1;
    end
    if depth==depthmax, error('unexpected error, the object is not recognized'), end
    % first object position
    subplot(double(ax))
    switch o(1).type
        case 'text'
            oldunits = o(1).units;
            o(1).units = 'normalized';
            position = o(1).position;
            outerposition = NaN;
            o(1).units = 'inches';
            pos = [o(1).extent(1)-margin*o(1).extent(3) o(1).extent(2)-margin*o(1).extent(4) (1+2*margin)*o(1).extent(3)  (1+2*margin)*o(1).extent(4)];
            createaxes = true;
            o(1).units = 'data';
            yax = o(1).extent(2)+[0-margin 1+margin]*o(1).extent(4);
            xax = o(1).extent(1)+[0-margin 1+margin]*o(1).extent(3);
            axlim = [xax yax];
            o(1).units = oldunits; % expected to be 'normalized'
            scale = getprop(double(ax),axesscaleproperties);
        case 'axes'
            oldunits = o(1).units;
            o(1).units = 'normalized';
            position = o(1).position;
            outerposition = o(1).outerposition;
            o(1).units = 'inches';
            pos = o(1).outerposition; % % absolute position of the outer region
            createaxes = false;
            axlim = axis;
            o(1).units = oldunits; % expected to be 'normalized'
            scale = getprop(double(o(1)),axesscaleproperties);
        otherwise
            oldunits = ax.units;
            ax.units = 'normalized';
            position = ax.position;
            outerposition = ax.outerposition;
            ax.units = 'inches';
            pos = ax.position; % absolute position of the inner region
            createaxes = true;
            axlim = axis;
            ax.units = oldunits;
            scale = getprop(double(ax),axesscaleproperties);
    end
    colormapcopy = get(fig,'colormap');
    subplot(double(ax)), [az,el] = view;
    % copy objects
    newfig = figure('visible','off','colormap',colormapcopy); % a new hidden figure
    delete(get(newfig,'children')) % empty figure
    if createaxes % any object but not axes
        newax = axes('position',[0 0 1 1]); % position
        newo = copyobj(o,newax);
        if strcmp(o(1).type,'text'), set(newo,'units','data'), end
        if length(axlim)==4
            set(newax,'xlim',axlim(1:2),'ylim',axlim(3:4),'visible','off');
        else
            set(newax,'xlim',axlim(1:2),'ylim',axlim(3:4),'zlim',axlim(5:6),'visible','off');
            view([az el])
        end
        if saveasbitmap==2
            set(newax,'visible','off')
            hlabels=get(newax,labels); % handles
            delete([hlabels{:}])
        end
    else % it is an axes object
        newo = copyobj(o,newfig);
        if saveasbitmap==2
            set(newo,'visible','off')
            hlabels=get(newo,labels); % handles
            delete([hlabels{:}])
        else
            hlabels=get(o,labels); % handles
            poslabels = get([hlabels{:}],'position'); % positions
            set(newo,'units','normalized','outerposition',[0 0 1 1])
            for ilabels=1:length(labels), set(get(newo,labels{ilabels}),'position',poslabels{ilabels}); end
        end
    end
    set(newfig,'paperunits','inches','paperposition',[0 0 pos(3:4)],'colormap',colormapcopy)

    % print and save
    if saveasbitmap
        if printoptions, print_jpg(res,filenametmp,printoptions)
        else print_jpg(res,filenametmp), end
    else
        saveas(newfig,filenametmp,'fig')
    end
end

% load image, crop it and store it as a MAT file (for efficiency: fast and reliable)
if saveasbitmap
    im = imread(filenametmp);
    if removewhitespace || saveasbitmap==2 % crop required
        if numel(im)<1e6 % small image
            asx = sum(sum(255-im,3),1);
            asy = sum(sum(255-im,3),2);
        else % too large image: switch to an iterative algorithm
            resim = size(im);
            asx = zeros(1,resim(2)); for i=1:resim(2), asx(i) = sum(sum(255-im(:,i,:),3),1); end
            asy = zeros(resim(1),1); for i=1:resim(1), asy(i) = sum(sum(255-im(i,:,:),3),2); end
        end
        thresh = size(im)/10;
        ix = [find(asx>thresh(1),1,'first') find(asx>thresh(1),1,'last')]; % crop assuming a white background
        iy = [find(asy>thresh(2),1,'first') find(asy>thresh(2),1,'last')];
        im = im(iy(1):iy(2),ix(1):ix(2),:);
    end
    sizeim = size(im);
    delete(filenametmp) % purge (not required since stored in im)
else
    sizeim = NaN;
    im = filenametmp;
end
date = datestr(now);
generator = 'obj2im';
save(filename,'im','sizeim','resolution','position','outerposition','scale','date','generator','-V7') % Matlab 7 format (more efficient)

% Purge
if any(newfig), delete(newfig), end

% output
if nargout, filestr = filename; end