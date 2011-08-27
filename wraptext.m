function s=wraptext(s0,width,newline,tol,trim)
%WRAPTEXT return a wrapped text with arbitrary newline sequence and accepting tolerance
%   syntax: s=wraptext(s0 [,width,newline,tol,trim])
%       width: required width (default=32)
%     newline: newline sequence (default='<br />')
%         tol: +/-tolerance (default=5)
%        trim: flag (default=true), remove leading and trailing spaces

% Migration 2.0 - 08/05/2011 - INRA\Olivier Vitrac - rev.

if nargin<2, width = 32; end
if nargin<3, newline = '<br />'; end
if nargin<4, tol =5; end
if nargin<5, trim =true; end

if iscell(s0), s = cellfun(@(si) wraptext(si,width,newline,tol,trim), s0,'UniformOutput',false); return, end
n = length(s0);
if n<=width, s=s0; return, end
if trim, s=strtrim(s0); else s=s0; end
p = regexp(s,'[\s-_\<\>\(\)\{\}\[\]\:\,\;\.]');
[dmin,imin] = min(abs(p-width));
if dmin<=tol, p=p(imin); else p=width; end
s = [s(1:p) newline wraptext(s(p+1:end),width,newline,tol) ];
