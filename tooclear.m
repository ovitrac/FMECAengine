function darkercolor = tooclear(color,threshold,correction)
%TOOCLEAR fixes too clear nx3 colors
%   Syntax: darkercolor = tooclear(color [,threshold,correction])
%       color = nx3 rgb matrix
%       threshold = gray level threshold for too clear (default = 0.6, for publications)
%       correction = gray level shift to make it darker (default = .15)
%
% INRA\MS 2.1 - 01/05/2016 - Olivier Vitrac - rev.

%% arg check
if nargin<1, error('one argument is required'), end
if ~isnumeric(color) || size(color,2)~=3, error('color must be a nx3 numeric array'), end
if nargin<2, threshold = []; end
if nargin<3, correction = []; end
if isempty(threshold), threshold = .6; end
if isempty(correction), correction = .15; end

darkercolor = color;
istooclear = mean(color,2)>threshold;
darkercolor(istooclear,:) = max(0,color(istooclear,:)-correction);