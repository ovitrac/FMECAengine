function dispf(varargin)
%DISPF fast wrapper of disp(sprintf(...))
% see help on SPRINTF
% see also FPRINTF (the main difference is that LF is used after disp)

% MS 2.1 - 16/03/08 - INRA\Olivier Vitrac rev. 19/08/18

% Revision history
% 29/12/12 updated help
% 19/08/18 implement varargin{1} as a cell

if iscell(varargin{1})
    for i=1:numel(varargin{1})
        dispf(varargin{1}{i})
    end
else
    disp(sprintf(varargin{:})) %#ok<DSPS>
end