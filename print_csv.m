function f=print_csv(filename,destination,varargin)
%PRINT_CSV print figure data into a CSV file
%     syntax: print_csv(filename [,folder,'property','value'])
%           filename: file name including path or not
%             folder: destination folder (use filename path if included)
%       list of poperty/default value
%               separator: ','
%           textdelimiter: '"'
%     note: f = print_csv(...) returns full filename
%
%   CSV file is assembled as 
%   +------------+------------+------------+------------+ ...
%   | "1:xlabel" | "1:ylabel" | "2:xlabel" | "2:ylabel" | ...
%   +------------+------------+------------+------------+ ...
%   | Xdata(1,1) | Ydata(1,1) | Xdata(1,2) | Ydata(1,2) | ...
%   +------------+------------+------------+------------+ ...
%   | Xdata(2,1) | Ydata(2,1) | Xdata(2,2) | Ydata(2,2) | ...
%   +------------+------------+------------+------------+ ...
%   |    ...
%   where Xdata and Ydata include the line data of all subplots
%
% See also: print_pdf, print_png, print_tif, gcfd, tab2csv

%INRA\MS 2.1 - 05/04/2015 - Olivier Vitrac - rev. 19/09/15

% Revision history
% 06/09/15 accept labels spread over multiple lines
% 19/09/15 use legends in titles (see how gcfd stores legend text in UserData)

% default
ext = '.csv';
default = struct(...
               'separator', ',',...
               'textdelimiter', '"',...
               'formattitle','%s ::%s' ...
               );

%% arg check
if nargin<1, error('one argument is required'), end
if nargin<2, destination = ''; end
[p,n] = fileparts(filename);
if isempty(destination)
    if isempty(p), destination = pwd; else destination = p; end
end
if ~exist(destination,'dir'), error('the destination folder ''%s'' does not exist',destination), end
filename = fullfile(destination,[n ext]);
options = argcheck(varargin,default);

%% extraction of data
data = gcfd;
if isempty(data) || ~isstruct(data) || ~isfield(data,'line'), error('no line data to export in CSV'), end
ndata = numel(data); blocksize = zeros(1,ndata); for i=1:ndata, blocksize(i) = length(data(i).line); end
idata = expandmat(1:ndata,blocksize);
xlabels = {data(idata).Xlabel}; inolabel = find(cellfun(@isempty,xlabels)); xlabels(inolabel) = arrayfun(@(i) sprintf('%d:xlabel',i),inolabel,'UniformOutput',false);
ylabels = {data(idata).Ylabel}; inolabel = find(cellfun(@isempty,ylabels)); ylabels(inolabel) = arrayfun(@(i) sprintf('%d:ylabel',i),inolabel,'UniformOutput',false);
lines = [data.line]; nlines = length(lines); nrows = arrayfun(@(i) length(lines(i).Xdata),1:nlines); m = max(nrows);

%% make table (cell array with labels as column title)
table = repmat({''},m+1,nlines);
for i=1:nlines
    icol = 1 + 2*(i-1);
    if isfield(lines,'UserData'), prefix = lines(i).UserData; else prefix=''; end
    if isempty(prefix) && nlines>4, prefix = sprintf('%02d',i); end
    if ~isempty(prefix)
        table{1,icol}   = sprintf(options.formattitle,printtxt(prefix),printtxt(xlabels{i}));
        table{1,icol+1} = sprintf(options.formattitle,printtxt(prefix),printtxt(ylabels{i}));

    else
        table{1,icol} = printtxt(xlabels{i});
        table{1,icol+1} = printtxt(ylabels{i});
    end
    table(2:nrows(i)+1,icol) = num2cell(lines(i).Xdata(:),2);
    table(2:nrows(i)+1,icol+1) = num2cell(lines(i).Ydata(:),2);
end
 
%% export table
tab2csv(table,filename,options)
fileinfo(filename)
if nargout, f = filename; end

%% prinvate function
function s=printtxt(s0)
if iscellstr(s0)
    s = regexprep(sprintf('%s\n',s0{:}),'\n','');
else
    s = s0;
end