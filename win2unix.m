function okout = win2unix(filename)
%WIN2UNIX converts WINDOWS ASCII file (with LF+CR) into UNIX ASCII file (with LF)
%
% Rem:
% the original file is renamed with a prefix 'MS_"
% the modified file has the same name as the original one.
%
% example:
% for i=1:6,win2unix(sprintf('C:\\Program Files\\Accelrys\\MS Modeling 3.0\\Gateway\\root_default\\dsd\\jobs\\Temp\\J%d.car',i)), end

% MS-MATLAB 1.0 - 09/01/04 - Olivier Vitrac - rev. 3/11/14

% Revision history
% 03/02/10 make compatible for linux
% 13/02/10 fix +h for windows
% 03/11/14 fix fileparts (version removed)

% Definitions
pref = 'MS_';
preftmp = 'tmp_';
ok = 1;

% file check
[pathstr,name,ext] = fileparts(filename);
destination = fullfile(pathstr,[pref name ext]);
destinationtmp = fullfile(pathstr,[preftmp name ext]);
if isempty(pathstr), pathstr = cd; end
if ~exist(filename,'file'),error('''%s'' does not exist in ''%s''',[name ext],pathstr), end
if exist(fullfile(pathstr,[pref name ext]),'file')
    dispf('''%s'' has already been converted',name)
    return
end

% load
data = textread(filename,'%s','delimiter','\n');

% conversion
fid = fopen(destinationtmp,'w');
if fid<=0,
    ok = 0;
else
    for i=1:length(data)
        fprintf(fid,'%s\n',data{i});
    end
    fclose(fid);
end

% rename
if ok
    if ~isunix, fileattrib(filename,'+h');  end                  %eval(sprintf('! attrib "%s" -H',filename))
    movefile(filename,fullfile(pathstr,[pref name ext]));       %eval(sprintf('! rename "%s" "%s"',filename,[pref name ext]))
    movefile(destinationtmp,fullfile(pathstr,[name ext]));      %eval(sprintf('! rename "%s" "%s"',destinationtmp,[name ext]))
    if ~exist(fullfile(pathstr,[pref name ext]),'file')
        error('unable to convert ''%s'' in ''%s''',[name ext],pathstr)
    end
end

% output
if nargout
    okout= ok;
elseif ok<=0
    disp('unable to convert')
end