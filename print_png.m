function print_png(resolution,fichier,chemin,options)
%PRINT_JPG  print active window as a PNG image
%   print_jpg(resolution,fichier[,chemin,options])
%   print_jpg(resolution,fichier,options)
%       options doit �tre pr�c�d� de '-' (ex. '-opengl')

% Woodox 1.0 - 28/01/01 - Olivier Vitrac - 13/08/07

% revision
% 25/07/07 add options and fix chemin option
% 13/08/07 fix path ambiguity when both fichier and chemin contain a path

if nargin<3, chemin = ''; end
if nargin<4, options = ''; end
if any(chemin), 
    if chemin(1)=='-', options = chemin; chemin = ''; end
    [pathstr,name,ext] = fileparts(fichier);
else
    [chemin,name,ext] = fileparts(fichier);
end
if isempty(chemin), chemin = cd; end
if ~exist(chemin,'dir'), error('the path ''%s'' does not exist',chemin), end
fichier = [name ext];
print (['-r' int2str(resolution)],'-dpng',options, fullfile(chemin,fichier))