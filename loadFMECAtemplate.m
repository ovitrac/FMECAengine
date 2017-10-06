function dbforsim = loadFMECAtemplate(filename,varargin)
% loadFMECAtemplate loads and tranform FMECAengine template created for 3SInPack into "sim" db
% excel file of template '...\INRA\Projects\3SinPack2016\template_FMECAengin\template_FMECAengine_V6.xlsm'
% sheetnames required:
%   'geometry'
%     'system'
%       'step'
% 'substances'
%'composition'
%          'K'
%        'CHI'
%
% EX. db = loadFMECAtemplate('C:\Data\Olivier\INRA\Projects\3SinPack2016\template_FMECAengine\template_FMECAengine_V6.xlsm');
% 
% INRA\Mai Nguyen, Olivier Vitrac
% 15/03/2016

% default definition
sheetnames = {'geometry' 'system' 'step' 'substances' 'composition' 'K' 'CHI'}; nsheets = length(sheetnames);
hearderbysheet = [8 1;     8  1;    6 1;      6 1;        7 1;     7 1;  7 1];
quantitytype = struct('Thickness', 'length',...
                  'Density', 'density',...
                  'Grammage', 'grammage',...
                  'Time','time',...
                  'Temperature','Ctemperature',...
                  'SML', 'weightconcentration',...
                  'Concentration', 'weightconcentration'); % correspondance beetween quantity in template and quantity in FMECAunit
density_default = struct(... SI units
                    'LDPE'                   ,[910 92],...
                    'LLDPE'                  ,[910 925],...
                    'HDPE'                   ,[942 965],...
                    'PP'                     ,[850 910],...
                    'PPrubber'               ,[850 910],...
                    'OPP'                    ,[850 910],...
                    'PS'                     ,[1040 1060],...
                    'HIPS'                   ,[1040 1060],...
                    'SBS'                    ,[1040 1060],...
                    'PET_sup_Tg'             ,[1370 1400],...
                    'PET_inf_Tg'             ,[1370 1400],...
                    'PBT'                    ,[1310 1340],...
                    'PEN'                    , 1000,...
                    'PA6'                    , 1130,...
                    'PA6_6'                  , 1140,...
                    'PA12'                   , 1020,...
                    'PVCrigid'               ,[1350 1420],...
                    'pPVC'                   ,[1350 1700],...
                    'AdhesiveNaturalRubber'  ,1000 ,...
                    'AdhesiveSyntheticRubber',1000,...
                    'AdhesiveEVA'            ,1000,...
                    'AdhesiveVAE'            ,1000,...
                    'AdhesivePVAC'           ,1000,...
                    'AdhesiveAcrylate'       ,1000,...
                    'AdhesivePU'             ,1000,...
                    'Paper'                  ,500 ,...
                    'Cardboard_polarmigrant' ,400,... polar substances
                    'Cardboard_apolarmigrant',400 ... hydrophobic substances
                        );
% arg check
if nargin<1, error('one argument is at least required: please give a fullname of template file'), end
if ~ischar(filename), error('the first argument must be a string'); end

% load file
% db.system = xlstblread(fullfile(local,templatefilename),'system',6,'headerrowindex',6);
for i = 1:nsheets
    db.(sheetnames{i}) = loadodsprefetch(filename,'sheetname',sheetnames(i),'headers',hearderbysheet(i,:),'prefetchsheet',1);
end
alllayerindex = db.system.Index(2:end);
alllayercode = db.system.Code(2:end);
allsubstancecode = db.substances.Code;
for s = sheetnames
    fieldbysheet.(s{1}) = fieldnames(db.(s{1}))';
end
% restreint data
% system sheet
ilayerok = cellfun(@isnan,db.system.Type,'UniformOutput',false);
ilayerok = find(~arrayfun(@(m) unique(ilayerok{m}),(1:length(ilayerok))'));
for f = fieldbysheet.system
    db.system.(f{1}) = db.system.(f{1})(ilayerok);
end
nlayer = max(db.system.Index(ilayerok)); % without food
% step sheet
istepok = cellfun(@isnan,db.step.Name,'UniformOutput',false);
istepok = find(~arrayfun(@(m) unique(istepok{m}),(1:length(istepok))'));
for f = fieldbysheet.step
    db.step.(f{1}) = db.step.(f{1})(istepok);
end
nsteps = max(db.step.Index(istepok));
% substance sheet
isolok =  cellfun(@isnan,db.substances.CASnumber,'UniformOutput',false);
isolok = find(~arrayfun(@(m) unique(isolok{m}),(1:length(isolok))'));
for f = fieldbysheet.substances
    db.substances.(f{1}) = db.substances.(f{1})(isolok);
end
nsubstance = max(db.substances.Index(isolok)); 
% composition sheet
for f = fieldbysheet.composition
    db.composition.(f{1}) = db.composition.(f{1})(isolok);
end
fcompositionremoved = arrayfun(@(m) sprintf('Concentration_L%d',m),alllayerindex(~ismember(alllayerindex,db.system.Index(2:end))),'UniformOutput',false);
db.composition = rmfield(db.composition,fcompositionremoved);
fieldbysheet.composition = fieldnames(db.composition)';
% K sheet
for f = fieldbysheet.K
    db.K.(f{1}) = db.K.(f{1})(isolok);
end
fKremoved = arrayfun(@(m) sprintf('KFP%d',m),alllayerindex(~ismember(alllayerindex,db.system.Index(2:end))),'UniformOutput',false);
db.K = rmfield(db.K,fKremoved);
fieldbysheet.K = fieldnames(db.K)';
% CHI sheet
for f = fieldbysheet.CHI
    db.CHI.(f{1}) = db.CHI.(f{1})(isolok);
end
fchiremoved = arrayfun(@(m) sprintf('chiP%d',m),alllayerindex(~ismember(alllayerindex,db.system.Index(2:end))),'UniformOutput',false);
db.CHI = rmfield(db.CHI,fchiremoved);
fieldbysheet.CHI = fieldnames(db.CHI)';

% define type of each column to be checked
for s = sheetnames(2:end)
    for f = fieldbysheet.(s{1})
        if strcmpi(f,'Index') % numeric
            typetobechecked.(s{1}).(f{1}) = 'noNaNnumeric';
        elseif strcmpi(f,'CASnumber') % CAS
            typetobechecked.(s{1}).(f{1}) = 'noNanCAS';
        elseif any(~cellfun(@isempty, (cellfun(@(m) strcmpi(regexp(f{1},'^Concentration','match'),m), fieldnames(quantitytype), 'UniformOutput',false))))
            typetobechecked.(s{1}).(f{1}) = 'noNaNquantity';
        elseif any(strcmpi(f,fieldnames(quantitytype))) % quantity class
            if strcmpi(f,'Thickness') || strcmpi(f,'Time') || strcmpi(f,'Temperature')
                typetobechecked.(s{1}).(f{1}) = 'noNaNquantity';
            else
                typetobechecked.(s{1}).(f{1}) = 'quantity';
            end
        elseif any(regexp(f{1},'^K')) % K
            typetobechecked.(s{1}).(f{1}) = 'K';
        elseif any(regexp(f{1},'^chi')) % chi
            typetobechecked.(s{1}).(f{1}) = 'chi';
        else % String
            if strcmpi(f,'Name') || strcmpi(f,'Type')
                 typetobechecked.(s{1}).(f{1}) = 'noNaNstring';
            else
                typetobechecked.(s{1}).(f{1}) = 'string';
            end
        end
    end
end

% check data in template
for s = sheetnames(2:end)
    for f = fieldbysheet.(s{1})
        tmp = struct([]);
        for i = 1:length(db.(s{1}).Index)        
            if any(regexp(typetobechecked.(s{1}).(f{1}),'quantity')) % quantity
                if any(regexp(f{1},'Concentration')) % if conc
                    if any(regexp(typetobechecked.(s{1}).(f{1}),'noNaN'))
                        type = ['noNaN' quantitytype.Concentration];
                    else type = quantitytype.Concentration; 
                    end
                else
                    if any(regexp(typetobechecked.(s{1}).(f{1}),'noNaN'))
                        type = ['noNaN' quantitytype.(f{1})];
                    else type = quantitytype.(f{1});
                    end
                end
                try 
                    [tmp(1).isok{i,1}, tmp(1).value{i,1}] = parserconvert(type,db.(s{1}).(f{1}){i});
                catch ME
                    if strcmp(ME.message,'MISSING DATA')
                        error('MISSING DATA at table "%s", column "%s" and row "%d"',s{1},f{1},i);
                    else
                        error('"%s" is not accepted for "%s"\ncheck table "%s", column "%s" and row "%d"',db.(s{1}).(f{1}){i},f{1},s{1},f{1},i)
                    end                    
                end
            else
                try
                    [tmp(1).isok{i,1}, tmp(1).value{i,1}] = parserconvert(typetobechecked.(s{1}).(f{1}),db.(s{1}).(f{1})(i));
                catch
                    error('MISSING DATA at table "%s", column "%s" and row "%d"',s{1},f{1},i);
                end
            end
        end
        parser.(s{1}).(f{1}) = struct2structtab(tmp);
    end
end

% fill 'empty' data in template (iok = 0 and value = '') by valid inputs with inheritance/ equation/ expertise/ ...

% interpret D and K, chi values
% K values
Kvalues = zeros(nsubstance,nlayer,nsteps);
layertype = db.system.Type(2:nlayer+1);
for isub = 1:nsubstance
    for ilayer = 1:nlayer
        for istep = 1:nsteps
            model = parser.K.(sprintf('KFP%d',ilayer))(isub).value;
            if isnumeric(model)
                Kvalues(isub,ilayer,istep) = model;
            else
                Kvalues(isub,ilayer,istep) = FMECAKtemplate(layertype(ilayer),model,db.step.Temperature{istep},db.substances.CASnumber{isub});
            end
        end
    end
end
% chi values
CHIvalues = zeros(nsubstance,nlayer+1,nsteps);
layertype = db.system.Type(1:nlayer+1);
for isub = 1:nsubstance
    for ilayer = 1:nlayer+1
        for istep = 1:nsteps
            if ilayer==1
                model = parser.CHI.chiF(isub).value;
            else
                model = parser.CHI.(sprintf('chiP%d',ilayer-1))(isub).value;
            end
            CHIvalues(isub,ilayer,istep) = FMECACHItemplate(layertype(ilayer),model,db.step.Temperature{istep});
        end
    end
end
Kvaluesfromchi = zeros(nsubstance,nlayer,nsteps);
layertype = db.system.Type(2:nlayer+1);
for isub = 1:nsubstance
    for ilayer = 1:nlayer
        for istep = 1:nsteps
            model = parser.CHI.(sprintf('chiP%d',ilayer))(isub).value;
            if isnumeric(model)
                Kvaluesfromchi(isub,ilayer,istep) = CHIvalues(isub,ilayer+1,istep)/CHIvalues(isub,1,istep);
            else
                Kvaluesfromchi(isub,ilayer,istep) = FMECAKtemplate(layertype(ilayer),'',db.step.Temperature{istep},db.substances.CASnumber{isub},CHIvalues(isub,1,istep));
            end
        end
    end
end
% replace empty Kvalues by Kvaluesfromchi
izero = Kvalues==0;
Kvalues(izero) = Kvaluesfromchi(izero);
% D model
Dvalues = zeros(nsubstance,nlayer,nsteps);
layertype = db.system.Type(2:nlayer+1);
diffusionmodel = db.system.Diffusionmodel(2:nlayer+1);
for isub = 1:nsubstance
    for ilayer = 1:nlayer
        for istep = 1:nsteps
            Dvalues(isub,ilayer,istep) = FMECADtemplate(layertype(ilayer),diffusionmodel(ilayer),parser.step.Temperature(istep).value,db.substances.CASnumber(isub));
        end
    end
end

% build data entry for FMECAengine
fiedtoberemoved = {'Index','Code','Name','CASnumber'};
CPforsim = rmfield(parser.composition,fiedtoberemoved);
k = 0;
for isub = 1:nsubstance
    if isub == 1
       dbforsim = fmecaengine('fmecamainfile',{'nlayers',nlayer,'nsteps',nsteps,...
            'constructor','lFm',parser.system.Thickness(1).value,...
            'l',{[parser.system.Thickness(2:end).value]},...
            'D',arrayfun(@(m) Dvalues(isub,:,m),1:nsteps,'UniformOutput',false),...
            'KFP',arrayfun(@(m) Kvalues(isub,:,m),1:nsteps,'UniformOutput',false),...
            'CP',{cellfun(@(f) parser.composition.(f)(isub).value,fieldnames(CPforsim))'}});
    else
       dbforsim(length(dbforsim)+1:length(dbforsim)+nsteps) = fmecaengine('fmecamainfile',{'nlayers',nlayer,'nsteps',nsteps,...
            'constructor','lFm',parser.system.Thickness(1).value,...
            'l',{[parser.system.Thickness(2:end).value]},...
            'D',arrayfun(@(m) Dvalues(isub,:,m),1:nsteps,'UniformOutput',false),...
            'KFP',arrayfun(@(m) Kvalues(isub,:,m),1:nsteps,'UniformOutput',false),...
            'CP',{cellfun(@(f) parser.composition.(f)(isub).value,fieldnames(CPforsim))'}});
    end    
    for istep = 1:nsteps
        k = k + 1;
        dbforsim(k).idusercode = sprintf('%s_%s',parser.substances.Code(isub).value,parser.step.Name(istep).value);
        if istep > 1
            dbforsim(k).parent = dbforsim(k-1).idusercode; 
        else
            dbforsim(k).parent = NaN;
        end
        dbforsim(k).inherit = NaN;
        dbforsim(k)
        dbforsim(k).ts = parser.step.Time(istep).value;
        dbforsim(k).SMLkgm3 = parser.substances.SML(isub).value;
        dbforsim(k).Binounit = parser.step.Biotnumber(istep).value; %{cellfun(@(f) parser.K.(f)(isub).value,fieldnames(Kforsim))'},...
       
    end
end

end % end of main function

%%%%========== PRIVATE FUNCTIONS ==========================
function [li,value] = parserconvert(type,strtest,testonly)
% function to check the validity of data and convert by FMECAunit if needed
% NaN value (empty cell in table): if required data -> li = 0, value = '' and error
%            else             -> li = 0; value = '' 
% other values: valid data -> li = 1; value = data;
%               else       -> li = 0; value = data;

    if nargin < 2, error('2 arguments are required: please define the type and the string to be tested'); end
    if nargin < 3, testonly = false; end
    if ~iscell(type), type = {type}; end
    if ~iscell(strtest), strtest = {strtest}; end
    % check NaN value/ data
    if isnan(strtest{1})
        if any(regexp(type{1},'noNaN')) % required data, NaN not accepted
            li = false; value = '';
            error('MISSING DATA')
        else li = false; value = '';
        end
    % check rules
    else % others value/ data
        if strcmpi(regexp(type{1},'numeric','match'),'numeric')  % numeric
            li = isnumeric(strtest{1}); value = strtest{1};
        elseif strcmpi(regexp(type{1},'string','match'),'string')  % char
            li = ischar(strtest{1}); value = strtest{1};
        elseif strcmpi(regexp(type{1},'CAS','match'),'CAS')  % CAS number
            li = checkCAS(strtest); value = strtest{1};
        elseif strcmpi(type{1},'K')  % partion coefficient
            if ischar(strtest{1}) && strcmpi(strtest{1},'KairP')
                li = true;  value = strtest{1};
            elseif isnumeric(strtest{1})
                li = true; value = strtest{1};
            else li = false; value = '';
            end
        elseif strcmpi(type{1},'chi')  % chi
            if ischar(strtest{1}) && strcmpi(strtest{1},'RT')
                li = true; value = strtest{1};
            elseif isnumeric(strtest{1})
                li = true; value = strtest{1};
            else li = false; value = '';
            end
        else  % quantity value-unit
            try
                tmp = FMECAunit(regexprep(type{1},'noNaN',''),strtest{1});
            catch
                li = false; value = '';
                error('"%s" is not accepted for "%s" quantity',strtest{1},type{1})
            end
            if isnumeric(tmp), li = true; value = tmp; end
            if ~testonly, li = true; value = tmp; end
        end
    end  
end % end parserconvert function

function Dmodel = FMECADtemplate(type,model,T,solute,solutepolarity)
%FMECADtemplate transform inputs in template into diffusion data which can be interpreted by FMECAengine
% INPUTS
% type: string, 'air','polymer', 'cardboard' (type of materials)
% model: string, 'Piringer' or 'Fuller' (diffusion model according to material type)
% T: array, temperature in °C
% solute: string, CAS number
% default
    T_default = 25; %celsius
    solute_default = 'anisole'; 
    % argcheck
    if nargin < 1, error('1 argument is required: Please define material type'), end
    if nargin < 2, model = ''; end
    if nargin < 3, T = [] ; end
    if nargin < 4, solute = ''; end
    if nargin < 5, solutepolarity = ''; end
    if isempty(model) 
        switch type{1}
            case 'air'
                model = 'Fuller';
            case 'food'
                model = 'Biot';
            otherwise
                model = 'Piringer'; % error('this material type (''%s'') has not been configured',type)
        end
    end
    if isempty(T), T = T_default; end
    if isempty(solute), solute = solute_default; end
    dbsolute = load_chemspider(solute);
    if isempty(dbsolute), error('the substance ''%s'' is not found',solute), end
    if isempty(solutepolarity)
        if ~isempty(dbsolute.Properties.LogP)
            if dbsolute.Properties.LogP.value > 0, solutepolarity = 'apolar';
            else solutepolarity = 'polar';
            end
        else solutepolarity = 'apolar';
        end
    end       
    % replace general type by specific type if required
    % 'dumb polymer', 'cardboard' : general type
    if strcmpi(type,'dumb polymer'), type = 'LDPE'; end
    if strcmpi(type,'cardboard')
        if strcmpi(solutepolarity,'polar'),  type = 'Cardboard_polarmigrant'; 
        elseif strcmpi(solutepolarity,'apolar'), type = 'Cardboard_apolarmigrant'; end
    end
    % build expression for calculation of D according to model
    switch model{1}
        case 'Fuller'
            formula = dbsolute.QuickMass.formula;
            Dmodel = Dfuller(formula,T);
        case 'Piringer'
            M = str2double(dbsolute.QuickMass.molecularmass);
            Dmodel = Dpiringer(type,M,T);
    %     case 'Biot'
    %         Dmodel = @(foodtexture) 
        otherwise
            error('this diffusion model (''%s'') has not been configured',model)
    end 
end % end FMECADtemplate

function Kmodel = FMECAKtemplate(type,model,T,solute,chi)
%FMECAKtemplate transform inputs in template into K values  which can be interpreted by FMECAengine
% INPUTS
% type: material/ phase type
% model: string or array, 'KairP' or numeric
% T: string, temperature interpreted by FMECAunit
% solute: string, CAS number
% chi: chi to calculate K
% default
    T_default = 25; %celsius
    solute_default = 'anisole'; 
    chi_default = 0;
    if nargin < 1, error('1 argument is required: Please define material type'), end
    if nargin < 2, model = ''; end
    if nargin < 3, T = [] ; end
    if nargin < 4, solute = ''; end
    if nargin < 5, chi = []; end
    if isempty(model)
        switch type{1}
            case 'air'
                model = 'KairP';
            otherwise
                error('this material type (''%s'') has not been configured',type);
        end
    end
    if isempty(T), T = T_default; end
    if isempty(solute), solute = solute_default; end
    if isempty(chi),  chi = chi_default; end
    switch model
        case 'KairP'
            Kmodel = 1/FMECAKairP(solute,T,chi);
        otherwise
            error('this partition model (''%s'') has not been configured',model)
    end 
end

function CHImodel = FMECACHItemplate(type,model,T)
%FMECACHItemplate transform inputs in template into K values  which can be interpreted by FMECAengine
% INPUTS
% type: material/ phase type
% model: string or array, 'KairP' or numeric
% T: string, temperature interpreted by FMECAunit
% default
    T_default = 25; %celsius
    if nargin < 1, error('1 argument is required: Please define material type'), end
    if nargin < 2, model = ''; end
    if nargin < 3, T = [] ; end
    if isempty(model)
        switch type{1}
            case 'air'
                model = 'RT';
            otherwise
                model = 0; % general value of CHI for condesed phase (food and polymer, etc.) % error('this material type (''%s'') has not been configured',type);
        end
    end
    if isempty(T), T = T_default; end
    switch model
        case 'RT'
            CHImodel = FMECAkair(T);
        otherwise
            CHImodel = FMECAgpolymer(model); % error('this CHI model (''%s'') has not been configured',model)
    end 
end