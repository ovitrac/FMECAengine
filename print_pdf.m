function print_pdf(resolution,fichier,chemin,options)
%PRINT_PDF print interactively the current figure as PDF document
%   print_PDF([resolution,fichier,chemin,options])
%   print_PDF generates a standard pdf based on figure properties
%   print_PDF(resolution,fichier,options)
%   print_PDF(resolution,fichier,chemin,options)
%
%   options 'nocheck' forces printing without any dialog
%
%   default options = '-loose'
%   for 3D without PLOTOBJ, options = '-opengl' is recommened (see PRINT for details)

% Woodox 1.0 - 13/08/07 - Olivier Vitrac - rev. 22/07/11

% revision
% 25/07/07 add options and fix chemin option
% 13/08/07 fix path ambiguity when both fichier and chemin contain a path
% 11/12/07 add fileinfo, preview and dialog
% 24/06/09 add nocheck
% 12/12/09 add evince &  for unix systems (in replacement of winopen)
% 20/07/11 fix chemin with 'nocheck'
% 21/07/11 new fix for 'nocheck'
% 22/07/11 fix extension with 'nocheck' 


% definitions
resolution_default = 600; % PDF
ext_default = '.pdf';
currentfig = gcf;
options_default = '-loose';
pmode = get(currentfig,'paperpositionmode');

% arg check
if nargin<1, resolution = []; end
if nargin<2, fichier = ''; end
if nargin<3, chemin = ''; end
if nargin<4, options = ''; end
if isempty(resolution), resolution = resolution_default; end
if isempty(fichier), [~,fichier] = fileparts(get(currentfig,'filename')); end
if isempty(fichier), [~,fichier] = fileparts(get(currentfig,'name')); end
if isempty(fichier), fichier = sprintf('Figure_%d',currentfig); end
if isempty(options), options = options_default; end
if any(chemin) 
    if chemin(1)=='-', options = chemin; chemin = ''; end
    [pathstr,name,ext] = fileparts(fichier);
    if ~isempty(pathstr) && exist(fullfile(chemin,pathstr),'dir')
        chemin = fullfile(chemin,pathstr);
    end
else
    [chemin,name,ext] = fileparts(fichier);
end
if isempty(chemin), chemin = pwd; end
if ~exist(chemin,'dir'), error('the path ''%s'' does not exist',chemin), end
if ~strcmp(ext,'.pdf'), ext = ext_default; end
fichier = fullfile(chemin,[name ext]);

% printing with no check
if strcmpi(options,'nocheck'), print(['-r' int2str(resolution)],'-dpdf','',fichier), return, end

% printing with controls
printon = true;
ok = ~exist(fichier,'file');
while ~ok
    answer = questdlg({sprintf('the file ''%s'' already exist',[name ext]),sprintf('in ''%s''',chemin)},'Overwrite an existing PDF file','overwrite','new file','cancel','overwrite');
    if strcmp(answer,'new file')
        [fichier,chemin] = uiputfile('*.pdf','Choose a new filename for your pdf');
        [chemin,name,ext] = fileparts(fullfile(chemin,fichier));
        if ~strcmp(ext,'.pdf'), ext = ext_default; end
        fichier = fullfile(chemin,[name ext]);
        ok = ~exist(fichier,'file');
    else
        printon = ~strcmp(answer,'cancel');
        ok = true;
    end
end
if printon && ~strcmp(pmode,'auto')
    answer=questdlg({ sprintf('Current ''paperpositionmode'' is ''%s''',pmode)
            ' '
            'If you set it as ''auto'', you can change the paper layout'
            'by resizing the figure and you can check the result by refreshing'
            'the preview in the printpreview panel.'
            'With ''auto'' printing starts only after you close the the preview panel'
            'You do not need to change the options in the preview panel.'
            ' '
            'Select an option' },'PAPERPOSITIONMODE: auto or not ?','auto','no change','cancel','no change');
    if strcmp(answer,'auto')
        set(currentfig,'paperpositionmode','auto');
    elseif strcmp(answer,'cancel')
        printon = false;
    else
        dispf('no change, paperpositionmode: %s',pmode)
    end
end
 
if printon
    if get(currentfig,'paperpositionmode')
        disp('Printing starts only after you close the the preview panel...')
        uiwait(printpreview)
    end
    start = clock;
    dispf('printing ''%s''....',fichier)
    print(['-r' int2str(resolution)],'-dpdf',options, fichier)
    dispf('... end in %0.3g s',etime(clock,start))
    fileinfo(fichier)
    if isunix
       system(['evince ' fichier ' &']);
    else
        winopen(fichier)
    end
else
    disp('PRINT_PDF canceled')
end