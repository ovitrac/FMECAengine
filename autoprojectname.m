function key = autoprojectname(keylength,tobefield,lowercase)
%AUTOPROJECTNAME generates a the name of a project (default length = 5)
% syntax: key = autoprojectname([keylength,tobefield,lowercase])
%   keylength: length of project name (default = 5)
%   tobefield: true to have a valid field name
%   lowercase: true to mix upper and lower cases
%         key: project name including only characters: 'A'..'Z' and '0'..'9'

% MS 2.1 - 02/02/08 - INRA\Olivier  - rev. 28/12/12

% Revision history
% 28/12/2012 add tobefield, lowercase, fix random seed

% default
keylength_default = 5;
tobefield_default = false;
lowercase_default = false;
rand('twister', sum(100*clock))

% arg check
if nargin<1, keylength = []; end
if nargin<2, tobefield = []; end
if nargin<3, lowercase = []; end
if isempty(keylength), keylength = keylength_default; end
if isempty(tobefield), tobefield = tobefield_default; end
if isempty(lowercase), lowercase = lowercase_default; end


% keygen
key = char(47+unidrnd(36,[1 keylength]));
ind = key>57;
key(ind) = key(ind)+7;

% modificators
if tobefield && key(1)<65
    key(1) = unidrnd(26,1,1) + 64;
end
if lowercase
    ind = find(key>58);
    u = rand(1,length(ind))>.5;
    key(ind(u)) = key(ind(u)) + 32;
end
