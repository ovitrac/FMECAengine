function chemin = find_path_toolbox(toolbox)
% FIND_PATH_TOOLBOX retourne le chemin d'installation d'une toolbox
%		ex. chemin = find_path_toolbox(toolbox)

% Woodox 1.0 - 03/03/01 - Olivier Vitrac - rev. 02/02/08

% revision history
% 02/02/08 add Unix compatibility

% arg check
if nargin<1, toolbox = ''; end
if isempty(toolbox), [root,toolbox] = lastdir(pwd); end
sep = pathsep; %if isunix, sep = ':'; else sep = ';'; end

% proceed matlabpath
matlabpath = path;
ind = findstr(lower(matlabpath), lower(toolbox));
if ~any(ind)
   chemin='';
else
   encadrement = nearest_index(ind(1),[0 findstr(matlabpath,sep) length(matlabpath)+1]);
   chemin = matlabpath(encadrement(1)+1:encadrement(2)-1);
end


function k = nearest_index(i,table)
k(1) = max(table(find(table-i<0)));
k(2) = min(table(find(table-i>=0)));