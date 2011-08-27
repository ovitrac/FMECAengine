function pngtruncateim(imfile,flipon,margin)
%PNGTRUNCATEIM  crop and rotate PNG already saved images
%   Syntax: pngtruncateim(imfile [,flipon,margin])

% Migration 2.0 - 24/05/11 - INRA\Olivier Vitrac - rev. 27/05/11

% Revision history
% 27/05/11 guess extension

if nargin<2, flipon = false; end
if nargin<3, margin = 100; end
im = imread(imfile); siz = size(im);
imb = min(im,[],3);
lim = zeros(2,2);
dimlist = 1:2;
for dim = dimlist
    lim(dim,1) = find(min(imb,[],dimlist(mod(dim,2)+1))<255,1,'first');
    lim(dim,2) = find(min(imb,[],dimlist(mod(dim,2)+1))<255,1,'last');
end
lim(:,1) = max(lim(:,1)-margin,1);
lim(:,2) = min(lim(:,2)+margin,siz(1:2)');
im = im(lim(1,1):lim(1,2),lim(2,1):lim(2,2),:);
if flipon, im = flipdim(permute(im,[2 1 3]),2); end
ext = uncell(regexp(imfile,'\.([^\.]+)$','tokens'));
if isempty(ext)
    imwrite(im,imfile,'png');
else
    imwrite(im,imfile,ext{1});
end