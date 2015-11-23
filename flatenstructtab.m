function sflat = flatenstructtab(sarray,fieldpattern,maxdupilicates,emptynum,emptytxt)
%FLATENSTRUCTTAB flatens a structure array by duplicating accordingly the fields
%   syntax: s = structtab2struct(sarray [,fieldpattern,maxdupilicates,emptynum,emptytxt])
%           sarray: structure array
%     fieldpattern: pattern to suffix duplicated fields (default='d%02d')
%   maxdupilicates: maximum number of accepted duplicates, larger number create an error (default=10)
%         emptynum: value for missing number (default = NaN)
%         emptytxt: value for missing text (default = 'none')
%   
%
%   example with dissimilar elements:
%{
    a.f1 = [1;2]; a.f2 = {{'a1';'a2';'a3'};'b'}; a.f3 = {(1:3)'; 2:3};
    b.f1 = [3;4;5]; b.f2 = {'e'; {'f';'g'}; {'h' 'i'}}; b.f3 = {3; (4:5)' ; (40:60)'};
    af = flatenstructtab(struct2structtab(a));
    bf = flatenstructtab(struct2structtab(b),[],100);
    abf = flatenstructtab([struct2structtab(a);struct2structtab(b)],[],100);
%}
%
%   See also STRUCT2STRUCTTAB, STRUCTTAB2STRUCT, SUBSTRUCTARRAY, TAB2CSV


% MS 2.1 - 23/11/15 - INRA\Olivier Vitrac - rev.

% Revision history
% 23/11/15 release candidate

%% Default inputs
fieldpattern_default = '_%02d';
maxdupilicates_default = 10;
emptynum_default = NaN;
emptytxt_default = 'none';

%% arg check
if nargin<1, error('one argument is required'), end
if ~isstruct(sarray), error('the first argument must be a structure'), end
if nargin<2, fieldpattern = ''; end
if nargin<3, maxdupilicates = []; end
if nargin<4, emptynum = emptynum_default; end
if nargin<5, emptytxt = emptytxt_default; end
if isempty(fieldpattern), fieldpattern = fieldpattern_default; end
if isempty(maxdupilicates), maxdupilicates = maxdupilicates_default; end

%% count the number of elements
n = numel(sarray);
fnames = fieldnames(sarray);
nfnames = length(fnames);
len = zeros(nfnames,1);
for ifname=1:nfnames
    for i=1:n
        counts = numel(sarray(i).(fnames{ifname}));
        len(ifname) = max(len(ifname),counts);
        if counts>maxdupilicates
            dispf('\tthe length (=%d) of sarray(%d).%s exceeds maxdupilicates=%d',counts,i,fnames{ifname},maxdupilicates)
        end
    end
end
if any(len>maxdupilicates)
    error('the length of one or several fields exceeds mexduplicates=%d',maxdupilicates)
end

%% prefetch the final structure
nnewfnames = sum(len);
newfnames = cell(nnewfnames,1);
expandedifnames = expandmat(1:nfnames,len);
counts = zeros(nfnames,1);
for ifname = 1:nnewfnames
    counts(expandedifnames(ifname)) = counts(expandedifnames(ifname))+1;
    if len(expandedifnames(ifname))>1
        newfnames{ifname} = sprintf(['%s' fieldpattern],fnames{expandedifnames(ifname)},counts(expandedifnames(ifname)));
    else
        newfnames(ifname) = fnames(expandedifnames(ifname));
    end
end
sflat = repmat(cell2struct(repmat({[]},nnewfnames,1),newfnames),size(sarray));

%% populate all fields
for i=1:n
    counts = zeros(nfnames,1);
    for ifname = 1:nnewfnames
        counts(expandedifnames(ifname)) = counts(expandedifnames(ifname))+1;
        if length(sarray(i).(fnames{expandedifnames(ifname)}))>=counts(expandedifnames(ifname))
            val = sarray(i).(fnames{expandedifnames(ifname)})(counts(expandedifnames(ifname)));
            if iscell(val), val = val{1}; end
            sflat(i).(newfnames{ifname}) = val;
        end
    end
end

%% clean empty values
for ifname = 1:nnewfnames
    val = {sflat.(newfnames{ifname})};
    isvoid = cellfun(@isempty,val);
    if any(isvoid)
        isnum = cellfun(@isnumeric,val);
        if all(isnum)
            [sflat(isvoid).(newfnames{ifname})] = deal(emptynum);
        else
            [sflat(isvoid).(newfnames{ifname})] = deal(emptytxt);
        end
    end
end