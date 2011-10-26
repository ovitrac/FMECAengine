function txt=formatsci(num,varargin)
%FORMATSCI print numbers as scientific Tex strings
%   Syntax: txt=formatsci(num [,param1,value1,param2,value2])
%   implemented properties
%       power10: anonymous function to calculate the power of 10, default value = @(x) floor(log10(x))
%    texpattern: valid formatting tex string with sprintf placeholders, default value =  %0.3g\cdot10^{%d}

% MS 2.1 - 25/10/11 - INRA\Olivier Vitrac - rev.

% arg check
param_default = struct(...
    'power10',@(x) floor(log10(x)),...
    'texpattern','%0.3g\cdot10^{%d}' ...
    );
param = argcheck(varargin,param_default);
param.texpattern = regexprep(param.texpattern,'\\','\\\');
texwhenremeq1 = uncell(regexp(param.texpattern,'^.*(10.*)','tokens'));

% power 10 and remainder
p10 = param.power10(num);
rem = num./10.^(p10);

% print collected values
% txt = sprintf(param.texpattern,[rem;p10]);
% txt = arrayfun(@(x,y) sprintf(param.texpattern,x,y),rem,p10,'UniformOutput',false);
n = numel(num);
txt = cell(1,n);
for i=1:n
    if abs(rem(i)-1)<sqrt(eps)/1e3;
        txt{i} = sprintf(texwhenremeq1{1},p10(i));
    else
        txt{i} = sprintf(param.texpattern,rem(i),p10(i));
    end
end
if n<2, txt = txt{1}; end