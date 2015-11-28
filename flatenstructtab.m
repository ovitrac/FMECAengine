function sflat = flatenstructtab(sarray,fieldpattern,maxdupilicates,emptynum,emptytxt,expandtextaslistflag)
%FLATENSTRUCTTAB flatens a structure array by duplicating accordingly the fields
%   syntax: s = structtab2struct(sarray [,fieldpattern,maxdupilicates,emptynum,emptytxt,expandtextaslistflag])
%           sarray: structure array
%     fieldpattern: pattern to suffix duplicated fields (default='d%02d')
%   maxdupilicates: maximum number of accepted duplicates, larger number create an error (default=10)
%         emptynum: value for missing number (default = NaN)
%         emptytxt: value for missing text (default = 'none')
% expandtextaslistflag: expand text separated by ';' or '|' as a list (default = true)
%   
%
% NOTE: In the examples below, the structures a, b and c are converted in structure arrays (structab)
%       before flatening. It is the original intent of this function. The function will work without
%       this conversion but will lead to different results (two calls of flatenstruct are required).
%
%       The new TABLE object type in MATLAB does not distinguish the difference between structure array
%       and accumulated structures(see the result of struct2table with or without prior conversion
%       with struct2structab).
% =======================================================================================
%   EXAMPLE WITH DISSIMILAR ELEMENTS (use at the end struct2table(af), to get this result)
% >> original table (use at the end struct2table(atest), to get this result)
%     f1        f2             f3          char 
%     __    __________    ____________    ______
% 
%     1     {3x1 cell}    [3x1 double]    'text'
%     2     'b'           [1x2 double]    'text'
%
% >> flatened table (use at the end struct2table(atestf), to get this result)
%     f1    f2_01    f2_02     f2_03     f3_01    f3_02    f3_03     char 
%     __    _____    ______    ______    _____    _____    _____    ______
% 
%     1     'a1'     'a2'      'a3'      1        2          3      'text'
%     2     'b'      'none'    'none'    2        3        NaN      'text'
%
%{
    a.f1 = [1;2]; a.f2 = {{'a1';'a2';'a3'};'b'}; a.f3 = {(1:3)'; 2:3};
    atest = struct2structtab(a); [atest.char] = deal('text');
    atestf = flatenstructtab(atest);
%}
% =======================================================================================
%
%   EXAMPLE 2 USING A VERY LARGE TABLE (note: maxdupilicates=100)
%{
    b.f1 = [3;4;5]; b.f2 = {'e'; {'f';'g'}; {'h' 'i'}}; b.f3 = {3; (4:5)' ; (40:60)'};
    bf = flatenstructtab(struct2structtab(b),[],100);
    abf = flatenstructtab([struct2structtab(a);struct2structtab(b)],[],100);
%}
% =======================================================================================
%
%   EXAMPLE 3 EXPAND TEXT AS LIST
% % >> original table (use at the end struct2table(c), to get this result)
%                 f1                                  f2                  
%     ___________________________    _____________________________________
% 
%     'A;B|C'    'AA;BB'    'CCC'    [1x3 double]    [1x2 double]    [100]
%
% >> flatened table (use at the end struct2table(cf), to get this result)
%     f1_01    f1_02     f1_03     f2_01    f2_02    f2_03
%     _____    ______    ______    _____    _____    _____
% 
%     'A'      'B'       'C'         1        2        3  
%     'AA'     'BB'      'none'     10       20      NaN  
%     'CCC'    'none'    'none'    100      NaN      NaN  
%{
    c.f1 = {'A;B|C' 'AA;BB' 'CCC'}; c.f2={1:3 10:10:20 100};
    cf = flatenstructtab(struct2structtab(c));
    abf = flatenstructtab([struct2structtab(a);struct2structtab(b)],[],100);
%}
% =======================================================================================
%
%   See also STRUCT2STRUCTTAB, STRUCTTAB2STRUCT, SUBSTRUCTARRAY, TAB2CSV


% MS 2.1 - 23/11/15 - INRA\Olivier Vitrac - rev. 25/11/2015

% Revision history
% 23/11/15 release candidate
% 25/11/15 fix char type, improved help and examples, add expandtextaslistflag

%% Default inputs
fieldpattern_default = '_%02d';
maxdupilicates_default = 10;
emptynum_default = NaN;
emptytxt_default = 'none';
expandtextaslistflag_default = true;

%% arg check
if nargin<1, error('one argument is required'), end
if ~isstruct(sarray), error('the first argument must be a structure'), end
if nargin<2, fieldpattern = ''; end
if nargin<3, maxdupilicates = []; end
if nargin<4, emptynum = emptynum_default; end
if nargin<5, emptytxt = emptytxt_default; end
if nargin<6, expandtextaslistflag = []; end
if isempty(fieldpattern), fieldpattern = fieldpattern_default; end
if isempty(maxdupilicates), maxdupilicates = maxdupilicates_default; end
if isempty(expandtextaslistflag), expandtextaslistflag = expandtextaslistflag_default; end

%% count the number of elements
n = numel(sarray);
fnames = fieldnames(sarray);
nfnames = length(fnames);
len = zeros(nfnames,1);
for ifname=1:nfnames
    for i=1:n
        if expandtextaslistflag && ischar(sarray(i).(fnames{ifname}))
            sarray(i).(fnames{ifname}) = expandtextaslist(sarray(i).(fnames{ifname}));
        end
        if ischar(sarray(i).(fnames{ifname}))
            counts = 1;
        else
            counts = numel(sarray(i).(fnames{ifname}));
        end
        len(ifname) = max(len(ifname),counts);
        if counts>maxdupilicates
            dispf('\tthe length (=%d) of sarray(%d).%s exceeds maxdupilicates=%d',counts,i,fnames{ifname},maxdupilicates)
        end
    end
end
if any(len>maxdupilicates)
    error('the length of one or several fields exceeds maxduplicates=%d',maxdupilicates)
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
            if ischar(sarray(i).(fnames{expandedifnames(ifname)}))
                val = sarray(i).(fnames{expandedifnames(ifname)});
            else
                val = sarray(i).(fnames{expandedifnames(ifname)})(counts(expandedifnames(ifname)));
            end
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