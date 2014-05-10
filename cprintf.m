function sf = cprintf(format,varargin)
% CPRINTF sprintf for cells (sprintf('format',cell1,cell2,text1,numeric1,...)
% Syntax: sf = cprintf(format,s1,s2,s3...)
%         s1,s2,s3.. are cell arrays, strinfs
%         sf is the formated string
% Alternative syntaxes
%         cprintf(format,s1,s2,s3...) alone behaves as dispf (without LF)
%         use cprintf({format formatadd1 formatadd2...},s1,s2,s3,...sn) to add special symbols after sn
%         NB: format is apply to all s1..sn
%             formatadd is accepted only when cprintf is used as a substitute of dispf
%             the main iitend is to enable the insertion of '\n' at the end of a long series of display
%             example: cprintf({'\tAUTOPREFETCH:\t[%20s]\n' '\n'},substructarray(whos,'name'))

% MS 2.1 - 24/12/12 - INRA\Olivier Vitrac - rev. 29/12/12

% Revision history
%  29/12/12 improved help, alternative syntaxes

if nargin<2, error('2 arguments are at least required'), end
if ~ischar(format) && ~iscellstr(format), error('format must be a string or cell array of strings'),end


% fix arguments: column vector
string = varargin;
for i=1:length(string)
    if ischar(string{i}), string{i} = cellstr(string{i}); end
    if iscell(string{i}), string{i} = string{i}(:); end
end
if any(cellfun(@iscell,string)), string = cat(1,string{:}); end

if nargout
    sf = sprintf(format,string{:});
else
    if iscell(format)
        fprintf(format{1},string{:});
        for i=2:length(format);        
            fprintf(format{i});
        end
    else
        fprintf(format,string{:});
    end
end