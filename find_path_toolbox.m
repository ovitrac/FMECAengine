function tbpath = find_path_toolbox(toolbox)
%FIND_PATH_TOOLBOX returns the path of a toolbox (rewritten version)
%   syntax: tbpath = find_path_toolbox(toolbox)
%   if no match is not found, a second search is performed while ignoring case
%   in case of multiple match the first is retrieved

% MS 2.1 - 03/03/01 - INRA/Olivier Vitrac - rev. 12/10/11

% Revision history
% 12/10/11: completely rewritten version with regular expressions to gain robustness, help in english

% arg check
if nargin<1 || nargin>1, error('one argument is required'), end
if ~ischar(toolbox), error('the first argument must be a string'); end
if isempty(toolbox), error('empty string'); end

% extract all paths and toolbox names
listtb_fullpath = regexp(path,pathsep,'split');
listtb = uncell(regexp(listtb_fullpath,['^.*\' filesep '(.*)$'],'tokens'));

% search the matching toolbox
found = find(~cellfun('isempty',regexp(listtb,['^' toolbox '$'])));
if isempty(found) % swicth to 'ignore case'
    found = find(~cellfun('isempty',regexpi(listtb,['^' toolbox '$'])));
end

% outputs
if any(found)
    tbpath = listtb_fullpath{found(1)};
    if length(found)>1
        warning('%s toolboxes match ''%s''',length(found),toolbox) %#ok<WNTAG>
    end
else
    tbpath = '';
end



%%%%%%%%%%%%%%%
% OLD CODE
%%%%%%%%%%%%%%%
% function chemin = find_path_toolbox(toolbox)
% % FIND_PATH_TOOLBOX retourne le chemin d'installation d'une toolbox
% %		ex. chemin = find_path_toolbox(toolbox)
% 
% % Woodox 1.0 - 03/03/01 - Olivier Vitrac - rev. 02/02/08
% 
% % revision history
% % 02/02/08 add Unix compatibility
% 
% % arg check
% if nargin<1, toolbox = ''; end
% if isempty(toolbox), [root,toolbox] = lastdir(pwd); end
% sep = pathsep; %if isunix, sep = ':'; else sep = ';'; end
% 
% % proceed matlabpath
% matlabpath = path;
% ind = findstr(lower(matlabpath), lower(toolbox));
% if ~any(ind)
%    chemin='';
% else
%    encadrement = nearest_index(ind(1),[0 findstr(matlabpath,sep) length(matlabpath)+1]);
%    chemin = matlabpath(encadrement(1)+1:encadrement(2)-1);
% end
% 
% 
% function k = nearest_index(i,table)
% k(1) = max(table(find(table-i<0)));
% k(2) = min(table(find(table-i>=0)));