function print_jpg(resolution,fichier,chemin,options)
%PRINT_JPG  crée un fichier jpeg99 à partir de la fenêtre active
%   print_jpg(resolution,fichier[,chemin,options])
%   print_jpg(resolution,fichier,options)
%       options doit être précédé de '-' (ex. '-opengl')

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
print (['-r' int2str(resolution)],'-djpeg99',options, fullfile(chemin,fichier))