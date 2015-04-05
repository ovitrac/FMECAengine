function txt=formatsci(num,varargin)
%FORMATSCI print numbers as scientific Tex strings
%   Syntax: txt=formatsci(num [,param1,value1,param2,value2])
%
%   Pair properties
%       power10: anonymous function to calculate the power of 10, default value = @(x) floor(log10(x))
%    texpattern: valid formatting tex string with sprintf placeholders, default value =  %0.3g\cdot10^{%d}
%                Note that '\' will be automatically protected
%                TIP: use mindigits when textpattern use %0.nf instead of %0.ng
%       pattern: alternative pattern to be compared when economoy is used
%     mindigits: minimum of digits to display (default = 0), best use with %0.nf
%
%   Keywords to force encony size representation
%       txt = formatsci(...,'economy') or formatsci(...,'eco')
%       
%   Example: formatsci([0.5e-15 12 1 1000 100000,-45,-89000,1e2],'eco') gives:
%               '5\cdot10^{-16}'    '12'    '1'    '10^{3}'    '10^{5}'    '-45'    '-8.9\cdot10^{4}'    '100'
%            without 'eco' or 'economy' the result is:
%               '5\cdot10^{-16}'    '1.2\cdot10^{1}'    '10^{0}'    '10^{3}'    '10^{5}'    '-4.5\cdot10^{1}'    '-8.9\cdot10^{4}'    '10^{2}'
%
%   Advanced example (for reports and publications)
%   The results are:
%         '10^{-3}'    '0.0314'    '1'    '1.41'    '3.14'    '10'    '31.4'    '100'    '500'
%         '10^{-3}'    '0.03'    '1.00'    '1.41'    '3.14'    '10^{1}'    '31.42'    '10^{2}'    '500.00'
%         '1.00\cdot10^{-3}'    '3.14\cdot10^{-2}'    '1.00'    '1.41'    '3.14'    '1.00\cdot10^{1}'    '3.14\cdot10^{1}'    '1.00\cdot10^{2}'    '5.00\cdot10^{2}'
%{
        numtotest = [1e-3 1e-2*pi 1 sqrt(2) pi 10 10*pi 100 500]
        formatsci(numtotest,'texpattern','%0.3g\cdot10^{%d}', 'pattern','%0.3g','eco') % shortest representation
        formatsci(numtotest,'texpattern','%0.2f\cdot10^{%d}', 'pattern','%0.2f','eco') % fixed digits instead
        formatsci(numtotest,'texpattern','%0.2f\cdot10^{%d}', 'pattern','%0.2f','eco','mindigits',2) % force a minimum of digits
%}

% MS 2.1 - 25/10/11 - INRA\Olivier Vitrac - rev. 17/03/15


% revision history
% 09/01/12 fix 0, add negative numbers
% 27/01/12 add pattern and economy
% 29/03/12 remove 10^0 fro results
% 08/12/13 improve output
% 14/03/15 force case
% 17/03/15 add 'eco' as a shorthand for 'economy)
% 17/03/15 major fix for formats including %0.nf instead of %0.nf (see advanced example)
% 17/03/15 fix formatsci(100,'eco') returns 100 instead of 10^2


% default values
param_default = struct(...
    'power10',@(x) floor(log10(x)),...
    'texpattern','%0.3g\cdot10^{%d}',...
    'pattern','%0.4g' ,...
    'mindigits',0 ...
    );
keywordlist = {'economy' 'eco'};

% arg checks
param = argcheck(varargin,param_default,keywordlist,'case');
param.texpattern = regexprep(param.texpattern,'\\','\\\');
texwhenremeq1 = uncell(regexp(param.texpattern,'^.*(10.*)','tokens'));

% power 10 and remainder
sgn = sign(num);
num = sgn.*num;
p10 = param.power10(num);
rem = num./10.^(p10);

% print collected values
% txt = sprintf(param.texpattern,[rem;p10]);
% txt = arrayfun(@(x,y) sprintf(param.texpattern,x,y),rem,p10,'UniformOutput',false);
n = numel(num);
txt = cell(1,n);
for i=1:n
    if num(i)==0
        txt{i} = '0';
    elseif abs(p10(i))~=0
        if (abs(rem(i)-1)<sqrt(eps)/1e3) && (param.mindigits==0);
            txt{i} = sprintf(texwhenremeq1{1},p10(i));
        else
            txt{i} = sprintf(param.texpattern,rem(i),p10(i));
        end
    else
        txt{i} = sprintf(param.pattern,rem(i));
    end
    if param.economy || param.eco
        tmp  = sprintf(param.pattern,num(i));
        tmpnum = tmp(tmp>='0' & tmp<='9');
        tmpnumsignificant = tmpnum(find(tmpnum>'0',1,'first'):end);
        tmpisallzero = isempty(tmpnumsignificant); %tmpisallzero = all(char(regexp(tmp,'\d','match'))=='0');
        tmptex = regexprep(txt{i},{'10\^','\^|{|}|\\cdot','^e'},{'e','','1e'});
        if (length(tmp)<=length(tmptex)) && ~any(tmp=='e')
            if ((~tmpisallzero) && length(tmpisallzero)>=param.mindigits) || num(i)==0
                txt{i} = tmp;
            end
        end
    end
    if sgn(i)<0, txt{i}=['-' txt{i}]; end
end
if n<2, txt = txt{1}; end