function dispf(varargin)
%DISPF wrapper of disp(sprintf(...))
% see help on SPRINTF

% MS 2.1 - 16/03/08 - INRA\Olivier Vitrac rev.

% Revision history

disp(sprintf(varargin{:}))