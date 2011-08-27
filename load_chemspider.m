function data = load_chemspider(mol,varargin)
%LOAD_CHEMSPIDER load CHEMSPIDER properties
%   syntax: data = load_chemspider(mol)
%           data = load_chemspider(mol,'property1',value1,'property2',value2,...)
%   INPUTS:
%      mol: string or nx1 cell array of strings containing molecules, CAS,...
%      recognized properties/values
%         destination: destination folder (default=tempdir)
%           thumbnail: flag (default=true), thumbnail image of the chemcial structure
%           structure: flag (default=true), 640x480 PNG image of the chemcial structure
%              follow: flag (default=false), true to force to follow all references found
%                csid: flag (default=false), true if mol contain CSID (to be used internally)
%   OUTPUTS:
%       data: nx1 structure array with fields
%            CSID: number, Chemispider identifier
%           InChI: string 
%        InChIKey: string
%          SMILES: string
%             CAS: string or cell array of strings
%        Synonyms: cell array of strings
%      Properties: structure with named fields (Averagemass,...VapourPressure) and containing a structure
%                  with fields: value (numeric), unit and name
%  UserProperties: structure array
%       Thumbnail: filename of the PNG thumbnail
%       structure: filename of the PNG stucture image
%             url: ChemSpider URL
%    urlstructure: ChemSpider URL for the structure
%
%
%   SEE ALSO: LOAD_NIST, LOAD_NIST_IR, LOAD_NCBI, LOAD_NCBISTRUCT, LOAD_CHEMINDUSTRY

% MS 2.1 - 21/05/11 - INRA\Olivier Vitrac - rev. 22/08/11

%Revision history
% 22/05/11 vectorization, minor bugs
% 13/06/11 set 'not found'
% 20/06/11 add iupac, remove tags <p>, <wbr>, <nobr>
% 13/07/11 fix requests based on CAS numbers starting with many 0
% 21/07/11 add QuickMass
% 22/08/11 update QuickMass according to modifications introduced by ChemSpider in August 2011

% default
keyword = {'thumbnail' 'structure'};
default = struct('thumbnail',true,'structure',true,'destination',tempdir,'csid',false,'follow',false);

% Configuration
token = '30c079d0-cedb-42a7-be40-6e8f9f2c0d75'; % Olivier Vitrac account
engine = 'http://www.chemspider.com/Chemical-Structure.%s.html';
imgengine = 'http://www.chemspider.com/ImagesHandler.ashx?id=%s&w=640&h=480';
tag2remove = {'<wbr\s*/?>' '</?nobr>' '</?p>' '</?span.*?>' '</?wbr>' '</?strong>' '</?sup>' '</?sub>' '\r|\n' '\s+' '\&\#(\d+)\;'};
tagreplacement = {'','','','','','','','','',' ','${char(str2double($1))}'};
num = '^\s*([+-]?\s*\d+\.?\d*[eEdD]?[+-]?\d*)'; % number

%arg check
if nargin<1, error('one argument is required'); end
options = argcheck(varargin,default,keyword,'property');
if ~exist(options.destination,'dir'), error('Destination folder ''%s'' does not exist',options.destination), end
if ~exist(fullfile(find_path_toolbox('MS'),'@Search'),'file')
    dispf('WARNING:\tCHEMSPIDER Search service is being installed, please wait...')
    currentpath = pwd; cd(find_path_toolbox('MS')) %#ok<*MCCD>
    createClassFromWsdl('http://www.chemspider.com/Search.asmx?WSDL')
    cd(currentpath)
    dispf('\t CHEMSPIDER Search service has been installed.\n\tFollow this: <a href="http://www.chemspider.com/Search.asmx">link</a> for details')
end

% recursion
if iscell(mol)
    nmol = numel(mol);
    data = struct([]);
    for i=1:nmol
        dispf('CHEMSPIDER iteration %d/%d: %s',i,nmol,mol{i})
        tmp = load_chemspider(mol{i},options);
        if isempty(data), data = tmp; else data(end+1:end+length(tmp))=tmp; end
    end
    return
end

%% Simple Search
tstart = clock;
obj = Search; % Constructor
% fix request based on CAS numbers and starting with many 0 
if all(checkCAS(mol))
    moltmp = regexp(mol,'[1-9]{1}[0-9]{1,6}-\d{2}-\d','match');
    if isempty(moltmp), error('unable to interpret %s as a valid CAS number',mol); end
    if length(moltmp)>1
        dispf('WARNING: several CAS numbers have been found, only the first is used:')
        cellfun(@(cas) dispf('\t''%s'' matches ''%s''',cas,mol),moltmp)
    end
    mol = moltmp{1};    
end
if ~options.csid
    screen = dispb('','LOAD_CHEMSPIDER\tsimple search service initiated for ''%s''',mol);
    try
        csid = SimpleSearch(obj,mol,token);
    catch errtyp
        if strcmp(errtyp.identifier,'MATLAB:unassignedOutputs'),
            dispb(screen,'WARNING LOAD_CHEMSPIDER: unable to find any valid CSID for %s',mol);
            data = struct([]); return
        else
            rethrow(errtyp)
        end
    end
    if isempty(csid), data = []; dispf('\t''%s'' not found',mol); return, end
    if iscell(csid.int)
        dispb(screen,'LOAD_CHEMSPIDER\t %d assessions found for ''%s''',length(csid.int),mol); screen='';
        if options.follow
            data = load_chemspider(csid.int,argcheck({'csid',true},options));
            return
        else csid = csid.int{1};
        end
    else    csid = csid.int;
    end
    screen = dispb(screen,'LOAD_CHEMSPIDER\t get information for ChemSpiderID=%s (''%s'')',csid,mol);
else
    screen = dispb('','LOAD_CHEMSPIDER\tuse the following ChemSpiderID ''%s''',mol);
    if isnumeric(mol), mol = num2str(mol); end
    if isempty(regexp(mol,'^\s*\d+\s*$','match','once')), error('invalid CSID ''%s''',mol), end
    csid = mol;
end
nfo = GetCompoundInfo(obj,csid,token); % first information

%% Extract thumbnail
if options.thumbnail
    screen = dispb(screen,'LOAD_CHEMSPIDER\t load thumbnail for ChemSpiderID=%s (''%s'')',csid,mol);
    thumbnailfile = fullfile(options.destination,sprintf('%s.thumb.png',csid));
    base64string=GetCompoundThumbnail(obj,csid,token);
    encoder = org.apache.commons.codec.binary.Base64;
    img = encoder.decode(uint8(base64string));
    fid = fopen(thumbnailfile,'w'); fwrite(fid,img,'int8'); fclose(fid);
else
    thumbnailfile = '';
end

%% Extract details
url = sprintf(engine,csid);
screen = dispb(screen,'LOAD_CHEMSPIDER\t connects to main URL %s',url);
details = urlread(url);
urlstructure = sprintf(imgengine,csid);
if options.structure
    screen = dispb(screen,'LOAD_CHEMSPIDER\t save 2D structure as image from URL %s',urlstructure);
    structurefile = fullfile(options.destination,sprintf('%s.png',csid));
    urlwrite(urlstructure,structurefile);
else
    structurefile = '';
end
dispb(screen,'LOAD_CHEMSPIDER\t extraction of ChemSpiderID=%s (''%s'') completed in %0.4g s',csid,mol,etime(clock,tstart));

%% Synonyms and CAS
iupac = strtrim(regexprep(uncell(regexp(details,'<div id="iupac".*?>(.*?)</div>','tokens')),tag2remove,tagreplacement));
syn = strtrim(regexprep(uncell(regexp(details,'<p class="syn".*?>(.*?)</p>','tokens')),tag2remove,tagreplacement));
cas = regexp(syn,'([1-9]{1}[0-9]{1,6}-\d{2}-\d)','tokens');
cas = uncell(cas(~cellfun('isempty',cas)));
if ~iscellstr(cas), cas = uncell([cas{:}]); end % additional uncell
cas = unique(cas);

%% Quick Properties (to be modified each time ChemSpider change its HTML code)
quickformula = regexprep(uncell(regexp(details,'<a.*class="emp-formula".*?>(.*?)</a>','tokens')),tag2remove,tagreplacement);
% Before August 2011 (not working now since at least Aug 22, 2011)
% quickmass = regexprep(uncell(regexp(details,'<li class="quick-mass">(.*?)</li>','tokens')),tag2remove,tagreplacement);
% quickmassprop    = regexprep(strtrim(regexprep(uncell(regexp(quickmass,'<th class="prop_title".*?>(.*?)</th>','tokens')),[tag2remove {'</?a.*?>' ':'}],[tagreplacement {'' ''}])),'\s','');
% quickmassvalue   = cellfun(@(x) str2double(x),strtrim(regexprep(uncell(regexp(quickmass,'<td class="prop_value.*?">([0-9\.]*)[\sDa]*?</td>','tokens')),tag2remove,tagreplacement)),'UniformOutput',false);
% quickmass = cell2struct([quickformula;quickmassvalue],[{'formula'};quickmassprop],1);
% After August 2011
quickmass = regexprep(uncell(regexp(details,'<li class="_quick-mass">(.*?)</li>','tokens')),tag2remove,tagreplacement);
quickmassprop = uncell(regexp(quickmass,'^\s*(.*):','tokens')); quickmassprop = regexprep(quickmassprop,'\s+','_');
quickmassvalue = uncell(regexp(quickmass,':\s*(.*) Da','tokens')); quickmassvalue = str2double(quickmassvalue{1});
quickmass = cell2struct([quickformula;quickmassvalue],[{'formula'};quickmassprop],1);

%% Properties
prop    = strtrim(regexprep(uncell(regexp(details,'<(td|th) class="prop_title".*?>(.*?)</(td|th)>','tokens')),[tag2remove {'</?a.*?>' ':'}],[tagreplacement {'' ''}]));
prop    = prop(:,2); % remove td|th
value   = strtrim(regexprep(uncell(regexp(details,'<td class="prop_value.*?">(.*?)</td>','tokens')),tag2remove,tagreplacement));
[valnum,stop]=regexp(value,num,'tokens','end');
valnum  = cellfun(@(x) str2double(x),uncell(valnum),'UniformOutput',false);
valunit = cellfun(@(x,k) strtrim(x(k+1:end)),value,stop,'UniformOutput',false);
fprop = regexprep(prop,{'ACD/' '#' '\(' '\s|\)|\.'},{'' 'Num' 'AT' ''}); [~,iprop] = unique(fprop); 
prop = cell2struct(cellfun(@(v,u,n) struct('value',v,'unit',u,'name',n),valnum(iprop),valunit(iprop),prop(iprop),'UniformOutput',false),fprop(iprop),1);

%% User properties
propuser = uncell(regexp(details,'<div.*?class="user_data_property_header_div".*?>(.*?)</div>','tokens'));
if ~isempty(propuser)
    propusername = strtrim(regexprep(uncell(regexp(propuser,'<span.*?>(.*)</span>','tokens')),':',''));
    propvalue = strtrim(regexprep(propuser,'<span.*?>.*?</span>',''));
    if ~isempty(propvalue)
        isnum = find(~cellfun('isempty',regexp(propvalue,[num '\s*$'],'once')));
        tmp = cellfun(@(x) str2double(x),propvalue(isnum),'UniformOutput',false);
        valid = cellfun(@(x) ~isnan(x),tmp);
        propvalue(isnum(valid)) = tmp(valid); 
    end
    propuser = struct('name',propusername,'value',propvalue);
else
    propuser = struct('name',{},'value',{});
end

%% Assembling
data = catstruct(nfo,struct('Name',iupac,'CAS',{cas'},'Synonyms',{syn'},'QuickMass',quickmass,'Properties',prop,'UserProperties',propuser,...
    'Thumbnail',thumbnailfile,'structure',structurefile,'url',url,'urlstructure',urlstructure));
if length(data.CAS)==1, data.CAS = data.CAS{1}; end
if length(data.Synonyms)==1, data.Synonyms = data.Synonyms{1}; end
data.CSID = str2double(data.CSID);
