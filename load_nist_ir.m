function sp = load_nist_ir(mol)
% LOAD_NIST_IR donwload IR spectra from nist
%   syntax: sp = load_nist(mol)
%    mol = string, array of string or object cretaed with  LOAD_CHEMINDUSTRY, LOAD_NCBI LOAD_NIST LOAD_CHEMINDUSTRY
%
% Examples: 
%   sp=load_nist_ir('128-37-0') -> only IR data are retrieved
%   sp=load_nist_ir(load_nist('128-37-0')) -> IR data are added to common data
%
% Other examples:
%  128-37-0 -> BHT
%  123-28-4 -> PS800
%  112-84-5 -> Erucamide

%   See also: LOAD_CHEMINDUSTRY, LOAD_NCBI LOAD_NIST LOAD_CHEMINDUSTRY

% ACTIALNE-MATLAB 2.0 - 15/04/09 - Olivier Vitrac - Guillaume Gillet - rev. 04/04/11

% Revision History
% 13/07/09 fix a bug when a single FTIR specctrum is proposed without listing page
% 04/04/11 update links to FTIR spectra

% definitions
root = 'http://webbook.nist.gov';
URL = sprintf('%s/cgi/cbook.cgi?Name=%%s&Units=SI',root);
tmpfile = fullfile(tempdir,'load_nist_ir.JCAMP-DX');

% arg check
if nargin<1, error('You may specify at least one molecule name or a CAS number'), end
if iscellstr(mol)
    if length(mol)>1
        sp = cell(size(mol));
        for i=1:numel(mol), sp{i} = load_nist_ir(mol{i}); end
        return
    else
        mol = mol{1};
    end
end
if isstruct(mol)
    if ~isfield(mol,'CAS') || ~isfield(mol,'name'), error('invalid object, must be created with LOAD_NIST, LOAD_NCBI, LOAD_CHEMINDUSTRY'), end
    sp = mol;
    for i=1:length(mol)
        if checkCAS(mol(i).CAS), sp.IR = load_nist_ir(mol(i).CAS); else sp.IR = load_nist_ir(mol(i).name); end
    end
    return
end
screen = '';
host = localname;

[iscas,validcas] = checkCAS(mol); if iscas, mol = validcas{1}; end
page = urlread(sprintf(URL,urlencoding(mol)));
links = regexpi(page,'<a href="/cgi/cbook.cgi\?ID=C\d*?&amp;Units=SI&amp;Mask=80#IR-Spec">.*?</a>','match');
title = regexprep(regexpi(page,'<title>.*?</title>','match'),'<.*?>','','ignorecase');
if ~isempty(links) && ~( strcmpi(title,'Name Not Found') || strcmpi(title,'Search Results') ) % bad names, links are available
    redirecturl = [root regexprep(links{1},'^.*?"|".*','')];
    screen = dispb(screen,sprintf('redirect to %s',regexprep(links{1},'<.*?>','')));
    page = urlread(redirecturl);
    if ~any(regexp(page,'Name Not Found')) && ~isempty(page)
        links = regexpi(page,'<a href="/cgi/cbook.cgi\?ID=C\d*?&amp;Units=SI&amp;Type=IR-SPEC&amp;Index=\d+#IR-SPEC">.*?</a>','match');
        if isempty(links), links = {''}; end
        for i=1:length(links)
            if ~isempty(links{i})
                redirecturl = [root regexprep(links{i},'^.*?"|".*','')];
                screen = dispb(screen,sprintf('redirect to %s',regexprep(links{i},'<.*?>','')));
                page = urlread(redirecturl);
            end
            % before 04/04/11
%             downloadlink = regexpi(page,'<a href="/cgi/cbook.cgi/[0-9-]*?\IR.jdx\?JCAMP=C\d*&amp;Index=\d*&amp;Type=IR">spectrum</a>','match');
            % after 04/04/11
            downloadlink = regexpi(page,'<a href="/cgi/cbook.cgi\?JCAMP=C\d*&amp;Index=\d*&amp;Type=IR">spectrum</a>','match');
            redirecturl = [root regexprep(downloadlink{1},'^.*?"|".*','')];
            screen = dispb(screen,sprintf('redirect to %s',downloadlink{1},'<.*?>',''));
            if exist(tmpfile,'file'), delete(tmpfile), end
            urlwrite(redirecturl,tmpfile);
            screen = dispb(screen,'\tread data...');
            out = jcampread(tmpfile);
            screen = '';
            for f=setdiff(fieldnames(out)',{'Blocks' 'Notes'}), sp(i).(f{1}) = out.(f{1}); end
            for f=fieldnames(out.Blocks)', sp(i).(f{1}) = out.Blocks.(f{1}); end
            prop = {}; val = {};
            for j=1:size(out.Notes,1)
                if (length(out.Notes{j,1})<20) || (j==1)
                    prop{end+1} = out.Notes{j,1};
                    val{end+1} =  out.Notes{j,2};
                else
                    val{end} =  [val{end} '; ' out.Notes{j,1}];
                    if ~isempty(out.Notes{i,2}), val{end} =  [val{end} '; ' out.Notes{j,2}]; end
                end
            end
            prop=strrep(prop,' ','_');
            for j=1:length(val)
                if ~isempty(regexp(val{j},'^(\s*-?\d+[\.]?\d*[Ee]?[-+]?\d*\s*)$','once')),  val{j}=str2double(val{j}); end
            end
            sp(i).nfo = cell2struct(val,prop,2);
            sp(i).url = redirecturl;
            sp(i).engine = mfilename;
            sp(i).date = datestr(now);
            sp(i).host = host;
        end
    else
        screen = dispb(screen,'%s: not found with NIST',mol);
        sp = [];
    end
else
    screen = dispb(screen,'%s: not found with NIST',mol);
    sp = [];
end