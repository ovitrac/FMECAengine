function data = load_chemspider(mol,varargin)
%LOAD_CHEMSPIDER load CHEMSPIDER properties
%   syntax: data = load_chemspider(mol)
%           data = load_chemspider(mol,'property1',value1,'property2',value2,...)
%   INPUTS:
%      mol: string or nx1 cell array of strings containing molecules, CAS,...
%      recognized properties/values
%         destination: destination folder (default=tempdir)
%           thumbnail: flag (default=true), thumbnail image of the chemcial structure
%           structure: flag (default=true), 640x480 PNG image of the chemical structure (see imsize)
%              follow: flag (default=false), true to force to follow all references found
%                csid: flag (default=false), true if mol contains a valid CSID (to be used internally)
%             nocache: flag (default=false), true = do not use the cache
%        refreshcache: flag (default=false), true = force the cache to be refreshed
%        reshashcache: flag (default=false), true = resfresh the hash corresponding to cached files (similar to clear functions)
%              imsize: 2x1 array coding for the image size of 2D structures (default = [640 480])
%            autocrop: flag (default=false)
%         orientation: 'horiz' or 'horizontal', 'vert' or vertical', 'none' (default)
%         cachefolder: fullfile(rootdir(which(mfilename)),'cache.ChemSpider')
%         autocorrect: it should be m x 2 cell array to be used as regexprep(name,autocorrect{1},autocorrect{2}) (default='')
%                      example {'^n-(.*ane)$','$1'} replaces 'n-pentane' in 'pentane'
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
%             EPI: structure coding for the different sections of EPI calculations
%       Thumbnail: filename of the PNG thumbnail
%       structure: filename of the PNG stucture image
%             url: ChemSpider URL
%    urlstructure: ChemSpider URL for the structure
%           links: list of links to other databases
%
%   Example: load_chemspider('2-nitrophenol','imsize',[1600,1200],'autocrop',true,'orientation','horiz')
%
%   NOTE: http://esis.jrc.ec.europa.eu/ provides additional ways to validate the quality of data stored in ChemSpider
%
%
%   SEE ALSO: LOAD_NIST, LOAD_NIST_IR, LOAD_NCBI, LOAD_NCBISTRUCT, LOAD_CHEMINDUSTRY

% MS 2.1 - 21/05/11 - INRA\Olivier Vitrac - rev. 24/04/18

%Revision history
% 22/05/11 vectorization, minor bugs
% 13/06/11 set 'not found'
% 20/06/11 add iupac, remove tags <p>, <wbr>, <nobr>
% 13/07/11 fix requests based on CAS numbers starting with many 0
% 21/07/11 add QuickMass
% 22/08/11 update QuickMass according to modifications introduced by ChemSpider in August 2011
% 11/01/12 update help with ESIS website
% 15/06/12 update quicklinks (including QuickMass), add imsize
% 19/09/13 add autocrop
% 05/10/13 works again to retrieve properties, extract links
% 06/10/13 add orientation
% 18/12/14 fix empty values (rare cases)
% 04/04/15 add EPI section data, add cache capabilities (major update)
% 07/04/15 fix with no cached files
% 27/05/15 TRANSIENT VERSION: first implementation of new chemspider rules
% 02/06/15 IMPROVED VERSION with almost full recovery of all ancient features, quickproperties
% 23/06/15 FIXES for EPI data and USER PROPERTIES
% 07/02/16 CONVERT user properties into numbers and extract units (to be completed, see recognizedunits)
% 08/02/16 add mmHg as recognized unit
% 02/04/16 trim mol input
% 19/05/16 add MAXNAMES due to inconsistencies in names aggregated by ChemSpider
%          example: Tridecane is also named "Dodecane", see: http://www.chemspider.com/Chemical-Structure.11882.html
%          currently MAXNAMES is set to 15 corresponding to last valid name for tridecane (93924-07-3)
% 24/11/16 fix unit properties with no space before (unit)
% 25/11/16 remove <stong> tags in user properties (to be checked in the future)
% 02/08/17 add autocorrect
% 16/04/18 fix CAS numbers with 0000 as prefix
% 24/04/18 fix ;

persistent CACHEDfiles CACHEglobal

% default
root = rootdir(which(mfilename));
keyword = {'thumbnail' 'structure'};
default = struct('thumbnail',true,...
                 'structure',true,...
                 'destination','',...
                 'cachefolder',fullfile(root,'cache.ChemSpider'),...
                 'csid',false,...
                 'follow',false,...
                 'nocache',false,...
                 'refreshcache',false,...
                 'reshashcache',false,...
                 'imsize',[640 480],...
                 'autocrop',false,...
                 'orientation','none',...
                 'autocorrect',[]);
MAXNAMES = 15;
makehash = @(id) lower([{id.CSID};id.CAS;{id.SMILES};{id.InChI};{id.InChIKey};id.names(1:min(length(id.names),MAXNAMES))]);

% Configuration
token = '30c079d0-cedb-42a7-be40-6e8f9f2c0d75'; % Olivier Vitrac account
%engine = 'http://www.chemspider.com/Chemical-Structure.%s.html'; % before 19/09/2013
rooturl = 'http://www.chemspider.com';
engine = sprintf('%s/RecordView.aspx?id=%%s',rooturl);
imgengine = sprintf('%s/ImagesHandler.ashx?id=%%s&w=%%d&h=%%d',rooturl);
tag2remove = {'</?span.*?>' '<wbr\s*/?>' '</?nobr>' '</?p>' '</?button.*?>' '</?wbr>' '</?strong>' '</?sup>' '</?sub>' '\r|\n' '\s+' '\&\#(\d+)\;'};
tagreplacement = {'','','','','','','','','','',' ','${char(str2double($1))}'};
num = '^\s*([+-]?\s*\d+\.?\d*[eEdD]?[+-]?\d*)'; % number
recognizedunits = {'°C' 'g/L' 'g/mL' 'mg kg-1' 'mmHg'}; % TO BE COMPLETED (OV: 08/02/16)

%arg check
if nargin<1, error('one argument is required'); end
options = argcheck(varargin,default,keyword,'property');
cachefolder = options.cachefolder;
pngcachefolder = fullfile(cachefolder,'thumbs');
if isempty(options.destination), options.destination = pngcachefolder; end
if ~exist(cachefolder,'dir'), mkdir(cachefolder), end
if ~exist(pngcachefolder,'dir'), mkdir(pngcachefolder), end
if ~exist(options.destination,'dir'), error('Destination folder ''%s'' does not exist',options.destination), end
if ~exist(fullfile(root,'@Search'),'file')
    dispf('WARNING:\tCHEMSPIDER Search service is being installed, please wait...')
    currentpath = pwd; cd(root) %#ok<*MCCD>
    createClassFromWsdl('http://www.chemspider.com/Search.asmx?WSDL')
    cd(currentpath)
    dispf('\t CHEMSPIDER Search service has been installed.\n\tFollow this: <a href="http://www.chemspider.com/Search.asmx">link</a> for details')
end
if ~isempty(options.autocorrect)
    if ~iscellstr(options.autocorrect) || size(options.autocorrect,2)~=2, error('autocorrect must be a mx2 cell of strings'), end
end

%cache index
if isempty(CACHEDfiles) ||  isempty(CACHEglobal) || options.reshashcache
    dispf('CHEMISPIDER: cache rehash')
    CACHEDfiles = explore('*.mat',cachefolder,0,'abbreviate'); nchachedfiles = length(CACHEDfiles);
    tmphash = cell(nchachedfiles,2);
    for i=1:nchachedfiles
        load(fullfile(CACHEDfiles(i).path,CACHEDfiles(i).file),'-mat','id')
        tmphash{i,1} = makehash(id); %#ok<NODEF>
        tmphash{i,2} = zeros(length(tmphash{i,1}),1,'uint16')+i;
    end
    CACHEglobal = struct('hash',{cat(1,tmphash{:,1})},'idx',{cat(1,tmphash{:,2})});
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

%% Use Cache if enabled
lowermol = lower(strtrim(mol));
if ~isempty(options.autocorrect), lowermol = regexprep(lowermol,options.autocorrect{:,1},options.autocorrect{:,2},'ignorecase'); end
[iscas,cascleaned]=checkCAS(lowermol); if iscas, lowermol = cascleaned; end % added 16/04/2018
if ~options.nocache && ~isempty(CACHEglobal.hash) && ismember(lowermol,CACHEglobal.hash)
    imol = CACHEglobal.idx(find(ismember(CACHEglobal.hash,lowermol),1,'first'));
    cachedfile = fullfile(CACHEDfiles(imol).path,CACHEDfiles(imol).file);
    if exist(cachedfile,'file')
        load(cachedfile,'data')
        dispf('CHEMSPIDER reuses cached data for ''%s'' (date=%s)',mol,CACHEDfiles(imol).date)
        return
    else
        dispf('WARNING: the cache of chemspider has been modified outside load_chemspider()')
        dispf('\tset use load_chemspider(....,''reshashcache'',true) or clear function to reshash the cache')
    end
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
        if strcmp(errtyp.identifier,'MATLAB:unassignedOutputs')
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
screen = dispb(screen,'LOAD_CHEMSPIDER\t connects to the main URL %s',url);
details = urlread(url); % next engine returns a link "<h2>Object moved to here.</h2>"

% links in page (parsing)
% current parser recognize: chemspider, wikipedia, pdb, google, ncbi, msds, jrc links...
linksinpage = uncell(regexp(details,'href="(.*?)"','tokens'));
isforeignlink = cellfun(@isempty,regexp(linksinpage,'^/'));
chemspiderlinks = unique(regexprep(...
      linksinpage(~isforeignlink ... keep URLs starting with /
      & cellfun(@isempty,regexp(linksinpage,'\.aspx$')) ... remove the following extensions
      & cellfun(@isempty,regexp(linksinpage,'\.ashx\?')) ...
      & cellfun(@isempty,regexp(linksinpage,'\.ico$')) ...
      & cellfun(@isempty,regexp(linksinpage,'\.css$')) ...
      & cellfun(@isempty,regexp(linksinpage,'\.pdf$')) ...
      & cellfun(@isempty,regexp(linksinpage,'^/$')) ...
      & cellfun(@isempty,regexp(linksinpage,'^/blog')) ... remove additional
      & cellfun(@isempty,regexp(linksinpage,'^/rss')) ...
      & cellfun(@isempty,regexp(linksinpage,'^/ChemSpiderOpenSearch')) ...
      ),'^/',sprintf('%s/',rooturl)));
isexternalink = isforeignlink ...
    & cellfun(@isempty,regexp(linksinpage,'^javascript')) ... remove javascript links
    & cellfun(@isempty,regexp(linksinpage,'^#')) ... remove internal links
    & cellfun(@isempty,regexp(linksinpage,'^http://oas.rsc.org/|^http://www.rsc.org/|^http://my.rsc.org')) ... % remove add links
    & cellfun(@isempty,regexp(linksinpage,'^\'''));
iswiki = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://en.wiki'));
ispdb    = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://www.rcsb.org/pdb/'));
isebi    = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://www.ebi.ac.uk'));
isgoogle = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://www.google.com'));
ismsds   = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://msds'));
isalfa   = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://www.alfa.com'));
isjrc    = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://esis.jrc.ec.europa.eu'));
isajpcell   = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://ajpcell'));
isncbi   = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://www.ncbi.nlm.nih.go'));
ismerck  = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://www.merckmillipore.com'));
isdoi    = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://dx.doi.org'));
isepa    = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://www.epa.gov'));
isacd    = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://www.acdlabs.com/'));
ischem    = isexternalink & ~cellfun(@isempty,regexp(linksinpage,'^http://www.chemicalize.org'));
parsedlinks = struct(...
   'chemspider', {chemspiderlinks},...
    'wikipedia', {linksinpage(iswiki)},...
       'google', {linksinpage(isgoogle)},...     
          'jrc', {linksinpage(isjrc)},...
         'ncbi', {linksinpage(isncbi)},...
          'epa', {linksinpage(isepa)},...         
          'pdb', {linksinpage(ispdb)},...
          'ebi', {linksinpage(isebi)},...
         'msds', {linksinpage(ismsds)},...
         'alfa', {linksinpage(isalfa)},...
         'cell', {linksinpage(isajpcell)},...
      'acdlabs', {linksinpage(isacd)},...
  'chemicalize', {linksinpage(ischem)},...          
        'merck', {linksinpage(ismerck)},...
          'doi', {linksinpage(isdoi)},...
       'others', {linksinpage(isexternalink & ~iswiki & ~isgoogle & ~isjrc & ~isncbi & ~isepa & ~ispdb & ~isebi & ~ismsds & ~isalfa & ~isajpcell & ~isacd & ~ischem & ~ismerck)} ...
       );
for link=fieldnames(parsedlinks)', if length(parsedlinks.(link{1}))<=1, parsedlinks.(link{1}) = char(parsedlinks.(link{1})); end, end 
      
% to be fixed / 19/09/2013 /// BUG fixed by Chemspider 05/10/2013
% linksinpage = uncell(regexp(details,'href="(.*?)"','tokens'));
% if length(redirectionlink)~=1, error('ChemSpider changed its policy, contact <a href="mailto:olivier.vitrac@agroparistech.fr">the developer</a>'), end
% details = urlread(redirectionlink{1});

% 2D structure
urlstructure = sprintf(imgengine,csid,options.imsize(1),options.imsize(2));
if options.structure
    screen = dispb(screen,'LOAD_CHEMSPIDER\t save 2D structure as image from URL %s',urlstructure);
    structurefile = fullfile(options.destination,sprintf('%s.png',csid));
    urlwrite(urlstructure,structurefile);
    if options.autocrop, pngtruncateim(structurefile,0,0,0,0); end
    if ~strcmpi(options.orientation,'none')
        imsize = imfinfo(structurefile);
        if (~isempty(regexpi(options.orientation,'^horiz')) && imsize.Width<imsize.Height), rotation = +90;
        elseif (~isempty(regexpi(options.orientation,'^vert')) && imsize.Width>imsize.Height), rotation = -90;
        else rotation =0;
        end
        if rotation, imwrite(imrotate(imread(structurefile),rotation,'nearest','loose'),structurefile); end
    end
else
    structurefile = '';
end
dispb(screen,'LOAD_CHEMSPIDER\t extraction of ChemSpiderID=%s (''%s'') completed in %0.4g s',csid,mol,etime(clock,tstart));

%% Synonyms and CAS
% iupac start to be discontinued from april 2015 (include French and German one)
% iupac = strtrim(regexprep(uncell(regexp(details,'<div id="iupac".*?>(.*?)</div>','tokens')),tag2remove,tagreplacement));
iupac = strtrim(regexprep(uncell(regexp(details,'<span class="prop_title">Systematic name</span>.*?<p>(.*?)</p>','tokens')),tag2remove,tagreplacement));
syn = strtrim(regexprep(uncell(regexp(details,'<div class="syn".*?>(.*?)</div>','tokens')),tag2remove,tagreplacement)); % before 27/05/15: <p></p> used instead of <div>
if isempty(iupac)
    acd = strtrim(regexprep(syn(~cellfun(@isempty,regexp(syn,'\[ACD/IUPAC Name\]')) & cellfun(@isempty,regexp(syn,'French')) & cellfun(@isempty,regexp(syn,'German'))),tag2remove,tagreplacement));
    iupac = strtrim(regexprep(acd,{'<.*>' '\[.*\]'},''));
end
cas = regexp(syn,'([1-9]{1}[0-9]{1,6}-\d{2}-\d)','tokens');
cas = uncell(cas(~cellfun('isempty',cas)));
if ~iscellstr(cas), cas = uncell([cas{:}]); end % additional uncell
cas = unique(cas);

%% Quick Properties (to be modified each time ChemSpider change its HTML code)
quickstr = regexprep(uncell(regexp(details,'<ul.*?class="struct-props".*?>(.*?)</ul>','tokens')),tag2remove(2:end),tagreplacement(2:end));
quickstr = [quickstr{:}];
quickpropengine = @(req)  strtrim(char(regexprep(uncell(regexp(quickstr,...
        sprintf('<span class="prop_title">%s</span>(.*?)</span>',req),'tokens')),[tag2remove {'\s*Copy\s*' '\s*Copied\s*'}],[tagreplacement,{'' ''}])));
quickpropenginenolink = @(req) strtrim(regexprep(quickpropengine(req),'</?a.*?>',''));
quickprop = struct(...
    'name',quickpropengine('Systematic name'),...
    'ChemSpiderID',regexprep(quickpropengine('ChemSpider ID'),'[^0-9]',''),...
    'MolecularFormula',quickpropengine('Molecular Formula'),...
    'AverageMass',regexprep(quickpropengine('Average mass'),'\s*Da.*$',''),...    
    'MonoIsotopicMass',regexprep(quickpropengine('Monoisotopic mass'),'\s*Da.*$',''),...
    'SMILES',quickpropengine('SMILES'),...
    'InChi',quickpropenginenolink('Std. InChi'),...
    'InChiKey',quickpropenginenolink('Std. InChIKey') ...
    );
% Before May 2015
% quickformula = regexprep(uncell(regexp(details,'<a.*class="emp-formula".*?>(.*?)</a>','tokens')),tag2remove,tagreplacement);
% Before August 2011 (not working now since at least Aug 22, 2011)
% quickmass = regexprep(uncell(regexp(details,'<li class="quick-mass">(.*?)</li>','tokens')),tag2remove,tagreplacement);
% quickmassprop    = regexprep(strtrim(regexprep(uncell(regexp(quickmass,'<th class="prop_title".*?>(.*?)</th>','tokens')),[tag2remove {'</?a.*?>' ':'}],[tagreplacement {'' ''}])),'\s','');
% quickmassvalue   = cellfun(@(x) str2double(x),strtrim(regexprep(uncell(regexp(quick{mass,'<td class="prop_value.*?">([0-9\.]*)[\sDa]*?</td>','tokens')),tag2remove,tagreplacement)),'UniformOutput',false);
% quickmass = cell2struct([quickformula;quickmassvalue],[{'formula'};quickmassprop],1);
% After August 2011
% quickmass = regexprep(uncell(regexp(details,'<li class="_quick-mass">(.*?)</li>','tokens')),tag2remove,tagreplacement);
% quickmassprop = uncell(regexp(quickmass,'^\s*(.*):','tokens')); quickmassprop = regexprep(quickmassprop,'\s+','_');
% quickmassvalue = uncell(regexp(quickmass,':\s*(.*) Da','tokens')); quickmassvalue = str2double(quickmassvalue{1});
% quickmass = cell2struct([quickformula;quickmassvalue],[{'formula'};quickmassprop],1);
% After June 14, 2012
% quicklinks not used anymore
%quicklinks = regexprep(uncell(regexp(details,'<div class="quick-links">(.*?)</div>','tokens')),tag2remove,tagreplacement);
averagemass = quickprop.AverageMass;   %uncell(regexp(quicklinks,'Average mass: (.+?) Da','tokens'));
monomass = quickprop.MonoIsotopicMass; %uncell(regexp(quicklinks,'Monoisotopic mass: (.+?) Da','tokens'));
quickformula = quickprop.MolecularFormula;
if isempty(quickformula), quickformula = ''; end
if isempty(averagemass), averagemass = ''; end
if isempty(monomass), monomass = ''; end
quickmass = cell2struct({quickformula;averagemass;monomass},{'formula' 'molecularmass' 'monoisotopicmass'},1);

%% Properties
prop    = strtrim(regexprep(uncell(regexp(details,'<(td|th) class="prop_title".*?>(.*?)</(td|th)>','tokens')),[tag2remove {'</?a.*?>' ':'}],[tagreplacement {'' ''}]));
if ~isempty(prop)
    prop    = prop(:,2); % remove td|th
    value   = strtrim(regexprep(uncell(regexp(details,'<td class="prop_value.*?">(.*?)</td>','tokens')),tag2remove,tagreplacement));
    [valnum,stop]=regexp(value,num,'tokens','end');
    valnum  = cellfun(@(x) str2double(x),uncell(valnum),'UniformOutput',false);
    valunit = cellfun(@(x,k) strtrim(x(k+1:end)),value,stop,'UniformOutput',false);
    fprop = regexprep(prop,{'ACD/' '#' '\(' '\s|\)|\.'},{'' 'Num' 'AT' ''}); [~,iprop] = unique(fprop);
    if ~isempty(valnum)
        prop = cell2struct(cellfun(@(v,u,n) struct('value',v,'unit',u,'name',n),valnum(iprop),valunit(iprop),prop(iprop),'UniformOutput',false),fprop(iprop),1);
    else
        dispf('\tWARNING: %d ChemSpider properties found but no values are attached',length(prop))
    end
else
    prop = '';
end

%% User properties (modified by OV: 07/02/2016)
propuser = uncell(regexp(details,'<span.*?class="user_data_property_name".*?>([^<>]{1,100}):</span></a></h2>\s*<table style="display:none">(.*?)</table>','tokens'));
% propuser = uncell(regexp(details,'<div.*?class="user_data_property_header_div".*?>(.*?)</div>','tokens'));
if ~isempty(propuser) && size(propuser,1)>1
    propusername = regexprep(propuser(:,1),{'\s' '[\(\)]'},{'',''});
    propvalue = uncell(regexp(propuser(:,2),'<td valign="top">(.*?)</td>','tokens'));
    if ~isempty(propvalue) && iscell(propvalue{1}) % fix 23/06/2015
        propvalue = cellfun(@(x) [x{:}]',propvalue,'UniformOutput',false);
    end
    % read numeric value
    npropval = length(propvalue);
    [propvaluenum,propvalunit] = deal(cell(npropval,1));
    for i=1:npropval
        if ~iscell(propvalue{i}), propvalue{i} = {propvalue{i}}; end %#ok<CCAT1> % added on Nov 25, 2016
        tmp = strtrim(regexprep(propvalue{i},{'<span.*?>.*?</span>' '<strong.*?>.*?</strong>' '<a.*?>.*?</a>' '\(.*\)' '\[.*\]' '(-?\d+\.?\d*)-(-?\d+\.?\d*)'},{'' '' '' '' '' '$1 $2'}));
        ntmp = length(tmp); propvaluenum{i} = zeros(1,ntmp); propvalunit{i} = cell(1,ntmp);
        for j=1:ntmp
            propvaluenum{i}(j) = mean(str2double(uncell(regexp(tmp{j},sprintf('(%s)',num(2:end)),'tokens'))));
            itmpunit = find(~cellfun(@isempty,regexp(propvalue{i}{j},recognizedunits,'match')),1,'first');
            if any(itmpunit), propvalunit{i}(j) = recognizedunits(itmpunit); else propvalunit{i}{j}=''; end
        end
        propvaluenum{i} = propvaluenum{i}(~isnan(propvaluenum{i}));
    end
%     propvalue = strtrim(regexprep(propuser,'<span.*?>.*?</span>',''));
%     if ~isempty(propvalue)
%         isnum = find(~cellfun('isempty',regexp(propvalue,[num '\s*$'],'once')));
%         tmp = cellfun(@(x) str2double(x),propvalue(isnum),'UniformOutput',false);
%         valid = cellfun(@(x) ~isnan(x),tmp);
%         propvalue(isnum(valid)) = tmp(valid);
%     end
    try
        propuser = cell2struct(num2cell(struct('value',propvaluenum,'unit',propvalunit,'raw',propvalue),2),regexprep(propusername,'[\s-\\/\(\)\[\]]*',''));
    catch
        propuser = struct('name',propusername,'value',propvaluenum,'unit',propvalunit,'raw',propvalue);
    end
else
    propuser = struct('name',{},'value',{});
end

%% EPI section data (added 04/04/15)
EPIclean = @(s) regexprep(strtrim(regexprep(s,{'\(.*?\)','\[.*?\]','/','-',','},{'','','_','_','_'})),'\s','');
EPIproppattern = {... generic parser (update it if new needs)
    '^(?<prop>[\w\s]*)\s?\((?<unit>.*)\)\s*[:=]\s*(?<num>[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)\s*\((?<nfo>.*?)\)'
    '^(?<prop>[\w\s]*)\s\((?<nfo>.*)\)\s*[:=]\s*(?<num>[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)\s*(?<unit>.*)'
    '^(?<prop>[\w\s]*)\s*[:=]\s*(?<num>[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)\s*(?<unit>.*)'
    '^(?<prop>[\w\s]*)\s*\[(?<nfo>.*)\]\s*[:=]\s*(?<num>[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?)\s*(?<unit>.*)'
    }; nEPIproppattern = length(EPIproppattern);
tmp = uncell(regexp(details,'<div class="tab-content" id="epiTab">(.*?)</div>','tokens'));
if ~isempty(tmp), tmp = uncell(regexp(tmp{1},'<pre.*?>(.*?)</pre>','tokens')); end
%tmp = uncell(regexp(details,'Predicted data is generated using the US Environmental Protection Agency’s EPISuite.*?<pre.*?>(.*?)</pre>','tokens'));
if isempty(tmp)
    EPIsections = struct([]);
    dispf('WARNING: no EPI data')
else
    EPIrawdata = strtrim(regexprep(tmp{1},'<.*?>.*</.*?>',''));
    EPItmp = regexp(EPIrawdata,'\n\n','split')';
    nsections = length(EPItmp);
    [EPItileofsection,EPIsectionname] = deal(cell(nsections,1));
    for i=1:nsections
        EPItmp{i} = strtrim(regexp(EPItmp{i},'\n','split'));
        EPItmp{i} = EPItmp{i}(~cellfun(@isempty,EPItmp{i}));
        EPItileofsection{i} = strtrim(regexprep(EPItmp{i}{1},':$',''));
        EPIsectionname{i} = EPIclean(EPItileofsection{i});
    end
    % remove duplicated EPI data (fix 23/06/2015)
    if ~isempty(findduplicates(EPIsectionname))
        dispf('WARNING: EPI data are duplicated')
        [EPIsectionname,u] = unique(EPIsectionname);
        EPItileofsection = EPItileofsection(u);
        EPItmp = EPItmp(u);
         nsections = length(EPItmp);
    end    
    EPIsections = cell2struct(repmat({struct('title','','raw','')},nsections,1),EPIsectionname);
    for i=1:nsections
        EPIsections.(EPIsectionname{i}) = struct('title',EPItileofsection{i},'raw',{EPItmp{i}(2:end)'});
        % read section
        for j=2:length(EPItmp{i})
            tmp = ''; iEPIproppattern = 0;
            while isempty(tmp) &&  iEPIproppattern<nEPIproppattern
                iEPIproppattern = iEPIproppattern+1;
                tmp = regexp(EPItmp{i}{j},EPIproppattern{iEPIproppattern},'names');
            end
            if ~isempty(tmp)
                if ~isfield(tmp,'nfo'), tmp.nfo='none'; end
                EPIsections.(EPIsectionname{i}).(EPIclean(tmp.prop)) = ...
                    struct('value',str2double(tmp.num),'unit',tmp.unit,'method',tmp.nfo);
            end
        end
    end
end
%% Assembling
% quick properties added on June 2, 2015
if isempty(cas), cas = {'NaN'}; end
data = catstruct(nfo,struct('Name',iupac,'CAS',{cas'},'Synonyms',{syn'},'QuickMass',quickmass,'quickproperties',quickprop,'Properties',prop,'UserProperties',propuser,...
    'EPI',EPIsections,'Thumbnail',thumbnailfile,'structure',structurefile,'url',url,'urlstructure',urlstructure,'links',parsedlinks));
if isfield(data,'CAS') && length(data.CAS)==1, data.CAS = data.CAS{1}; end
if isfield(data,'Synonyms') && length(data.Synonyms)==1, data.Synonyms = data.Synonyms{1}; end
data.CSID = str2double(data.CSID);
if ~isfield(data,'structure'), data.structure = structurefile; end
if ~isfield(data,'url'), data.url = url; end
if ~isfield(data,'urlstructure'), data.urlstructure = urlstructure; end

%% cache update if required
cachedfilename = sprintf('%d.mat',data.CSID);
cachedfile = fullfile(cachefolder,cachedfilename);
if options.refreshcache || ~exist(cachedfile,'file')
    if ~isfield(data,'CAS'), data.CAS = 'NaN'; end
    if ~iscell(data.CAS), idCAS = {data.CAS}; else idCAS = data.CAS; end
    id = struct('CSID',sprintf('%d',data.CSID),'CAS',{idCAS(:)},'SMILES',data.SMILES,'InChI',data.InChI,'InChIKey',data.InChIKey,...
        'names',{strtrim(regexprep(data.Synonyms(:),{'<a.*?>.*?</a>','\[.*?\]'},{'' ''}))});
    save(cachedfile,'id','data')
    tmphash = makehash(id);
    if isempty(CACHEglobal.idx), imol = 1; else imol = CACHEglobal.idx(end)+1; end
    if isempty(CACHEDfiles)
        CACHEDfiles = explore(cachedfilename,cachefolder,0,'abbreviate');
    else
        CACHEDfiles(end+1) = explore(cachedfilename,cachefolder,0,'abbreviate');
    end
    CACHEglobal = struct('hash',{[CACHEglobal.hash;tmphash]},'idx',[CACHEglobal.idx;zeros(length(tmphash),1,'uint16')+imol]);
    dispf('LOAD_CHEMISPIDER: updated cache'), fileinfo(cachedfile)
end
