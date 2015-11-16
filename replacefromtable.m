function xr = replacefromtable(x,table)
%REPLACEFROMTABLE replaces any character (or words) by an other from a table of substitution
% SYNTAX: xr = replacefromtable(x,table)
%   x = char or cell array of strings
%   table = substitution table (mx2 char or cell array, empty value = character to be deleted)
%   xr = substitued string (char or cell array as x)

% MIGRADIM 1.0 - 08/03/03 INRA\Olivier Vitrac rev. 05/06/09

% revision history
% 05/06/09 word implementation

% Input control
if ~iscell(table), ischaron_table =1; else ischaron_table =0; end
if ~iscell(x), x = cellstr(x); ischaron_x =1; else ischaron_x =0; end
[m,n] = size(table);
if n>m && m>1
	[m,n] = size(table');
	table = table';
end

% Replace
for u = 1:length(x)
    xr{u} = x{u};
    for i = 1:m
        if ischaron_table
            lookfor = table(i,1); % for char
            replace = table(i,2);
        else
            lookfor = table{i,1}; % for cellstr
            replace = table{i,2};
        end
        if length(lookfor)<2
            j = find(xr{u}==lookfor);
            if any(j)
                if any(replace)
                    xr{u}(j) = replace;
                else
                    xr{u}(j) = ''; %xr{u}(setdiff(1:length(x),j));
                end
            end
        else % words
            xr{u} = regexprep(xr{u},lookfor,replace);
        end
    end
end

if ischaron_x, xr = char(xr); end