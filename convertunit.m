function varargout = convertunit(varargin)
%UNIT converts between units, also checks dimensions
%   convertunit(str1,str2) displays the conversion factor between str1 and str2.
%   examples: convertunit('N','kg.m/s2') or convertunit('MW','J/s')
%   note the units can be compound terms with powers and prefixes
%
%   convertunit(str1,str2,num) converts the specified number of units
%
%   f        = convertunit(...) returns the conversion as an output
%   [f,resu] = convertunit(...) returns additional information about the units
%
%   convertunit('report') displays a report of the code:
%                  units      (e.g. meter)
%                  unit kinds (e.g. length)
%                  prefixes   (e.g. M for mega=10^6)
%                  options    (e.g. 'abstemp' for absolute temperatures)
%
%   other options: list, report, show, verb, debug, eng, abstemp, showfrac
%
%TYPICAL EXAMPLES (not that scalar values can be replaced by vectors or arrays)
%   t=convertunit('day','s',10) %converts 10 days into seconds
%   v=convertunit('m/s','km/hour',20) %converts 20 m/s into km/h
%   l=convertunit('um','m',30) % converts 30 um into m
%   C=convertunit('mg/L','kg/m^3',40) % converts 40 mg/L (almost ppm) into kg/m3
%   T=convertunit('C','K',50,'abstemp') % converts 50°C into K (absolute temperature, not temperature difference)
%   A=convertunit('cm^2','m^2',60) % converts 60 cm2 into m2
%   D=convertunit('cm^2/s','m^2/s',7e-9) % converts 7e-9 cm2/s into m2/s
%
%EXAMPLES
%   convertunit('kg','lb')       %converts from kilos to pounds
%   convertunit('atm','lbf/in2') %converts from standard atmospheres to psi.
%   convertunit('lbf/lb','m/s2') %acceleration due to gravity in m/s^2
%   convertunit('lbf-in','N.m')  %moment conversion
%   convertunit('cc','m3',3.2)   %convert 3.2 cc to m3
%   convertunit('N.m','kg*m2/s^2')
%   convertunit('W','kW.hr/day') %kiloWatt-hours per day
%   convertunit('rpm','Hz')      %Hz has dimensions rev/second, not just 1/second
%
%UNITS
%   convertunit('report') displays all units recognised by this function
%   convertunit('list')   displays just their names
%
%   these may be compound units, e.g. N.m or lb/s
%   terms must be separated by one of . * / or space
%
%   each term may have a power, e.g m^3 or m3. the ^ is optional.
%
%   each term may have a scaling prefix, e.g. kW for kiloWatt
%   powers are applied to prefixes: e.g. km3 is km x km x km
%
%   divisors are applied to all following terms, e.g.
%   kg/m*s2 is kg^1 m^-1 s^-2 : note how the power of "s" is negative
%
%CHECKS
%   UNIT checks that input/output units have identical dimensions
%   e.g. that both are velocities, or mass/time, etc.
%   this is done by reducing each unit to basic "kinds" of unit, e.g length
%   convertunit('report') will display all the kinds of unit
%
%CUSTOM UNITS
%   it is possible to define additional custom units as inputs
%   these new units may not use the name/symbol of an existing unit
%   custom units may reference built-in units, or earlier custom units
%   custom units may be combined into 1xn structure arrays
%   example:
%   Y.name = 'year'; Y.symb = 'yr'; Y.size = 365.242; Y.expr = 'day';
%   convertunit('year','second') %this will fail, because "year" is unknown
%   convertunit('year','second',Y) %this will work
%   example:
%   Y.name = 'year'; Y.symb = 'yr'; Y.size = 365.242; Y.expr = 'day';
%   Z.name = 'LightYear'; Z.symb = 'ly'; Z.size = 3E8; Z.expr = 'yr*m/s';
%   convertunit('LightYear','m',Z,Y); %this will fail, because Z references Y
%   convertunit('LightYear','m',Y,Z); %this will work
%   example:
%   Y.name = 'year'; Y.symb = 'yr'; Y.size = 365.242; Y.expr = 'day';
%   Z.name = 'LightYear'; Z.symb = 'ly'; Z.size = 3E8; Z.expr = 'yr*m/s';
%   N = [Y Z];
%   convertunit('LightYear','m',N); %this uses a set of custom units
%
%OUTPUTS
%   first output is the conversion factor from first to second input:
%   factor = convertunit('lb','kg'); %factor = 0.4536
%   if an amount is received, first output is conversion of that many:
%   amt = convertunit('lb','kg',[3 4 5]); gives amt = [1.3608 1.8144 2.2680];
%
%   second output is additional information
%   [~,resu] = convertunit('W','kJ/hr',2)
%   resu.conv = 3.60 conversion factor
%   resu.in.symb: 'W'          %given name at input
%   resu.in.expr: 'P1.'        %kind of unit: power^1
%   resu.in.dims: 'l2.m1.t-3.' %dimensions length^2 * mass / time^3
%   resu.in.size: 1            %input (W) divided by base unit (kgm2/s3)
%   resu.in.amt : 2            %number of units input (1 by default)
%   resu.out.symb: 'kJ/hr'     %symbol from input
%   resu.out.expr: 'E1.t-1.'   %kind of unit: energy per time
%   resu.out.dims: 'l2.m1.t-3.' %dimensions must match input above
%   resu.out.size: 0.2778      %input (kJ/hr) divided by base unit (kgm2/s3)
%   resu.out.amt: 7.200 %conversion factor multiplied by amount received
%
%   if custom units are received, the following fields will be added:
%   resu.customUnits.name
%   resu.customUnits.symb
%   resu.customUnits.expr
%   resu.customUnits.size
%
%DISPLAY
%   by default, numbers are displayed using the standard %g format
%   Numbers close to a multiple of pi will be shown as pi x %g
%   'eng' option display all numbers using engineering format
%   'frac' format displays all dimensions as fractions
%   examples:
%   convertunit('rev/s','rad/s')          will display 2pi
%   convertunit('rev/s','rad/s','eng')    will display 6.28E0
%
%VERBOSITY
%   there are three levels of verbosity
%   'show'  displays the input/output unit dimesions and conversion factor
%           on by default when the function has no output(s)
%           off by default when the function has output(s)
%   'verb'  above plus shows lookup from unit to dimensions of:
%           (1) all the diffeent kinds of unit (e.g. power, force, area)
%           (2) input/output units
%   'debug' above plus shows each step of each lookup
%
%TEMPERATURES
%   by default, temperature conversion is by scale, not absolute
%   'abstemp' option converts absolute temperatures
%   will throw an error if you try to use 'abstemp' with a unit that isn't
%   a pure temperature
%   examples:
%   convertunit('C','F',3)              returns 5.4 farenheit per 3 celcius
%   convertunit('C','F',3,'abstemp')    returns 37.4 farenheit = 3 celcius
%   convertunit('C/s','F/s','abstemp'); returns an error
%
%version: 2015.01.29 - modified INRA\Olivier Vitrac - 01/04/15, 04/04/15

%% SETUP
if ~nargin, convertunit('report'), return, end
opt = sub_setup_opt(nargout,varargin{:}); %parse options
opt = sub_setup_prefix(opt);      %list of prefixes, e.g. M = mega = 1E6
opt = sub_setup_constant(opt);    %constants, e.g. g=9.81m/s
opt = sub_setup_kind(opt);        %kinds of unit, e.g. mass.
opt = sub_setup_base(opt);        %all base units with scale, eg lb, kg, g

%% SHOW UNIT REPORT, STOP IF NO CONVERSION REQUESTED
if opt.report; sub_report(opt); end
if opt.list;   sub_list(opt); end
if isequal(opt.in,'none'); return; end

%% DIMENSION LOOKUP AND CONVERSION FACTOR
%find the units in and out, and their conversion factor from base
opt.in = sub_get_unit_data(opt,opt.in);
opt.ud = sub_get_unit_data(opt,opt.ud);
opt.conv = opt.in.size/opt.ud.size;

%check unit dimensions match
sub_check_units(opt);

%warn if converting single power of units without the abstemp option
[nin,pin] = sub_parse_units(opt.in.dims);
if isequal(nin,{'T'}) && isequal(pin,1) && ~opt.abstemp
    warning([mfilename ':abstemp'],'This will convert the scale of temperatures, not absolute temperatures. You may want the ''abstemp'' option.');
end

%convert absolute temperatures, if requested
if opt.abstemp; opt = sub_temp_convert(opt); end

%% ASSIGN OUTPUT
if nargout;
    resu.out = opt.ud;
    resu.in = opt.in;
    resu.conv = opt.conv;
    if opt.abstemp;
        resu.in.amt = opt.in.size;
        resu.out.amt = opt.ud.size;
    else
        resu.in.amt = opt.amt;
        resu.out.amt = opt.amt.*resu.conv;
    end
    if ~isempty(opt.newunit); resu.customUnits = opt.newunit; end
    varargout{1} = resu.out.amt;
    varargout{2} = resu;
end
if ~nargout || opt.show; sub_show_units(opt); end
end

%% SUBROUTINES : SETUP
function opt = sub_setup_opt(nout,varargin)

%extract options from input
opts = {'show','verb','debug','report','eng','abstemp','frac','list'};
for i=1:numel(opts)
    opt.(opts{i}) = false;
end
%by default, show is true if no outputs from unit
if ~nout; opt.show = true; end

opt.newunit = [];
keep = true(size(varargin));
for i=1:numel(varargin)
    %deal with new units (defined as structures)
    if isstruct(varargin{i});
        rqfields = {'name','symb','size','expr'};
        found = isfield(varargin{i},rqfields);
        if ~all(found)
            disp(rqfields(~found));
            error('new units must have the above fields');
        end
        
        if isempty(opt.newunit); opt.newunit = varargin{i};
        else                     opt.newunit = [opt.newunit varargin{i}];
        end
        keep(i) = false;
        continue;
    end
    
    if isnumeric(varargin{i}); continue; end %skip numbers
    switch varargin{i}
        case opts; opt.(varargin{i}) = true; keep(i) = false;
    end
end
opt.verb = opt.debug || opt.verb;
varargin = varargin(keep);
any_opt = numel(varargin)<numel(keep);

%handle rest of input which should be 0-3 arguments:
switch numel(varargin)
    case 0; %only allow this if any options are true
        if ~any_opt; error('need inputs: either units or options'); end
        opt.in = 'none';
        opt.ud = 'none';
        opt.amt = nan;
    case 1;
        error('need 2 inputs: units in and out');
    case 2;
        opt.in = varargin{1};
        opt.ud = varargin{2};
        opt.amt = 1;
    case 3;
        %handle number not last
        if isnumeric(varargin{1}); varargin = varargin([2 3 1]); end
        if isnumeric(varargin{2}); varargin = varargin([1 3 2]); end
        opt.in = varargin{1};
        opt.ud = varargin{2};
        opt.amt = varargin{3};
    otherwise
        error('can only accept 0-3 inputs, plus options');
end

%check inputs are correct form
if ~ischar(opt.in); error('first input must be a string'); end
if ~ischar(opt.ud); error('second input must be a string'); end
if ~isnumeric(opt.amt); error('third input must be a number or option'); end

if numel(size(opt.in))>2 || size(opt.in,1)>1 || any(size(opt.in)==0)
    error('first input must be a 1xn string, cannot be empty');
end

if numel(size(opt.ud))>2 || size(opt.ud,1)>1 || any(size(opt.ud)==0)
    error('second input must be a 1xn string, cannot be empty');
end
end

function opt = sub_setup_kind(opt)
%creates data on "kinds" of unit
%opt = sub_setup_kind(opt)

opt.kind = sub_data_kind();
opt.kind.dims = sub_simplify(opt.kind.expr);
opt.kind.size = ones(size(opt.kind.dims));

%the subset of kind that is dims
opt.dims = sub_subset(opt.kind,opt.kind.isdim);

    function kind = sub_data_kind()
        b = {...
            'length',       'l','l',    'meter',...
            'time',         't','t',    'second',...
            'temperature',  'T','T',    'Celcius',...
            'charge',       'c','c',    'coulomb',...
            'infosize',     's','s',    'bit',...
            'mass',         'm','m',    'kg',...
            'rev',          'rev','rev','rev',...
            'area',         'A','l2',   'meter^2',...
            'volume',       'V','l3',   'meter^3',...
            'velocity',     'v','l/t',  'meter/s',...
            'acceleration', 'a','v/t',  'meter/second^2',...
            'force',        'F','m.a',  'Newton',...
            'energy'        'E','F.l',  'Joule',...
            'power',        'P','F.v',  'Watt',...
            'frequency',    'f','rev/t','Hertz',...
            'current',      'I','c/t',  'Ampere',...
            };
        
        %removed, now defined from power&current
        %'resistance',   'R','R',    'R',...
        
        kind.name = b(1:4:end);
        kind.symb = b(2:4:end);
        kind.expr = b(3:4:end);
        kind.unit = b(4:4:end);
        kind.numberof = numel(kind.name);
        kind.isdim = false(size(kind.name));
        kind.isdim(1:7) = true;
        
        %name is a human-friendly name
        %symb is the one-letter symbol for this kind (must be unique)
        %expr is the expression for this, using symbols from "symb" field
        %unit all units are scaled relative to this, so e.g. all masses are
        %relative to kg (NOT grammes!)
        %you cannot reference this table directly when calling the
        %function, because "length/time" is not a unit, but "m/s" is.
    end

    function a = sub_subset(a,flag)
        for f=fieldnames(a)'
            v = a.(f{1});
            if numel(v)==numel(flag);
                v = v(flag);
                a.(f{1}) = v;
            end
        end
    end
end

function opt = sub_setup_prefix(opt)
%sets up data on prefixes
%opt = sub_setup_prefix(opt)

opt.prefix.symb = 'numcHkMGT';
opt.prefix.powr = [-9 -6 -3 -2 2 3 6 9 12];
end

function opt = sub_setup_constant(opt)
%sets up conversion constants
%opt = sub_setup_con(opt)

%standard gravity, as used by e.g. lbf->N conversion
opt.con.g = 9.80665;
opt.con.abszero = -273.15;

%opt.con.c = 299792458; %speed of light in m/s

opt.con.C2K = @(t) t-opt.con.abszero;
opt.con.C2R = @(t) (t-opt.con.abszero)*1.8;
opt.con.C2F = @(t) t*1.8 + 32;

opt.con.K2C = @(t) t+opt.con.abszero;
opt.con.R2C = @(t) t/1.8+opt.con.abszero;
opt.con.F2C = @(t) (t-32)/1.8;
end

function opt = sub_setup_base(opt)
%set up data on base units
%opt = sub_setup_base(opt)

opt.unit = sub_base_data();
opt.unit = sub_base_newunit(opt.unit,opt);
opt.unit.numberof = numel(opt.unit.name);
opt = sub_base_convert(opt);

    function opt = sub_base_convert(opt)
        %for any units that are references to other units, convert to references to
        %the base units
        any_changed = true;
        if opt.verb;
            disp('=======================================');
            disp('== all units : convert to dimensions ==');
            disp('=======================================');
        end
        while any_changed;
            any_changed = false;
            for i=1:opt.unit.numberof
                
                u = opt.unit.expr{i};
                if opt.verb;
                    disp('--------');
                    fprintf(1,'unit %s with expression %s\n',opt.unit.symb{i},opt.unit.expr{i});
                end
                
                [opt.unit.expr{i},opt.unit.size(i)] = sub_dimensions(opt,opt.unit.expr{i},opt.unit.size(i),'base',false,opt.debug);
                
                if opt.verb
                    fprintf(1,'is dimensions %s scale %f\n',opt.unit.expr{i},opt.unit.size(i));
                end
                any_changed = any_changed || ~isequal(u,opt.unit.expr{i});
                %if any_changed; keyboard; end
            end
            
        end
        if opt.verb;
            disp('===============================')
            disp('== finding base units : DONE ==')
            disp('===============================');
        end
    end

    function out = sub_base_data()
        b = {'meter',   'm',    1,      'l', ...
            'second',   's',    1,      't', ...
            'second',   'sec',  1,      't', ...
            'celsius',  'C',    1,      'T', ...
            'degC',     '°C'     1,      'T', ...
            'farenheit','F',    5/9,    'T', ...
            'kelvin',   'K',    1,      'T', ...
            'rankine',  'Ra',   5/9,    'T', ...
            'coulomb',  'c',    1,      'c', ...
            'amp',      'A',    1,      'I', ...
            'ohm',      'R',    1,      'P/I2', ...
            'volt',     'V',    1,      'P/I', ...
            'bit',      'bit',  1,      's', ...
            'gram',     'g',    1E-3,   'm', ...
            'rev',      'rev',  1,      'rev',...
            'minute',   'min',  60,     't',...
            'hour',     'hr',   3600,   't', ...
            'day',      'day',  24,     'hr',...
            'week',     'week',  7,     'day',...
            'month',    'month', 30.5,  'day',...
            'year',     'year',  365.25,'day',...
            'ton',      'T',    1E3,    'm', ...
            'newton',   'N',    1,      'F', ...
            'litre',    'L',    1E-3,   'V', ...
            'cubicCentimeter','cc',   1E-6,   'V', ...
            'usGallon', 'gal', 1/264.172,'V',...
            'joule',    'J',    1,      'E', ...
            'watt',     'W',    1,      'P', ...
            'herz',     'Hz',   1,      'f', ...
            'rpm',      'rpm',  1,      'rev/min',...
            'degree',   'deg',  1/360,  'rev',...
            'radian',   'rad',  1/2/pi, 'rev',...
            'pascal',   'Pa',   1,      'N/l2', ...
            'pound',    'lb',   0.453592,'m', ...
            'kiloforce', 'kgf', opt.con.g, 'N',...
            'poundforce','lbf', 0.453592*opt.con.g,'N', ...
            'foot',     'ft',   0.3048, 'l',...
            'inch',     'in',   0.0254, 'l',...
            'wattHour','Wh',   3600,   'J',...
            'mile',     'M',   1609.344,'l',...
            'nauticalMile','NM',1852   ,'l',...
            'stoke',    'St',   0.01,    'l2/sec',...
            'psi',      'psi',  6894.75729, 'Pa',...
            'stdatmo' ,'atm',  101325,'Pa',...
            'inchesMercury','inHg',3386,'Pa',...
            'mmMercury','mmHg',133.32239,'Pa',...
            'bar','bar',1E5,'Pa'
            };
        out.name = b(1:4:end);
        out.symb = b(2:4:end);
        out.size = [b{3:4:end}];
        out.expr = b(4:4:end);
    end

    function U = sub_base_newunit(U,opt)
        pat = 'custom unit "%s" bad format: field "%s" must be a %s\n';
        pat2 = 'custom unit "%s" clashes: field "%s" value "%s" already exists\n';
        for ii = 1:numel(opt.newunit)
            new = opt.newunit(ii);
            %check the new unit is compatible
            if ~ischar(new.name);       error(pat,'','name','string'); end
            if ~ischar(new.expr);       error(pat,new.name,'expr','string'); end
            if ~ischar(new.symb);       error(pat,new.name,'symb','string'); end
            if ~isnumeric(new.size);    error(pat,new.name,'size','number'); end
            
            %warn about anything being overwritten
            if ismember(new.name,U.name); error(pat2,new.name,'name',new.name); end
            if ismember(new.symb,U.symb); error(pat2,new.name,'symb',new.symb); end
            
           
            ud = sub_base_unit(opt,new.expr);
            U.name{end+1} = new.name;
            U.symb{end+1} = new.symb;
            U.size(end+1) = ud.size*new.size;
            U.expr{end+1} = ud.dims;
            opt.unit = U;
        end
        
        
        if (numel(opt.newunit) && opt.show) || opt.debug || opt.verb
            fprintf('Using %d custom units:\n',numel(opt.newunit));
            pat = 'name "%s" expression 1 %s = %g x %s\n';
            for ii=1:numel(opt.newunit)
                fprintf(pat,...
                    opt.newunit(ii).name,...
                    opt.newunit(ii).symb,...
                    opt.newunit(ii).size,...
                    opt.newunit(ii).expr);
            end
            fprintf('\n');
        end
        
    end

end

%% UNITS IN/OUT
function out = sub_get_unit_data(opt,str)
%parse arbitrary unit to get dimensions and scale
%opt = sub_get_unit_data(opt,str)

%parse separators to all be .
loc = regexp(str,'-[a-z|A-Z]','start');
str(loc) = '.';

%first part : lookup into base units
out = sub_base_unit(opt,str);

if opt.verb;
    disp(['== input unit "' str '" : START ==']);
    fprintf(1,'original expression %s\n',str);
    fprintf(1,'->\n');
end

%second part : lookup into dimensions (subset of kinds)
[out.dims,s] = sub_dimensions(opt,out.expr,1,'kind',true ,opt.debug);
out.size = out.size * s;

if opt.debug; fprintf(1,'=== finished with this unit\n'); end
if opt.verb;
    fprintf(1,'new expression\n');
    temp = out;
    temp.symb = out.dims;
    sub_show_dimensions(opt,temp)
    fprintf(1,'scale %f\n',out.size);
    disp(['=== input unit "' str '" : DONE ===']);
end
end

%% PARSING OF UNITS AND DIMENSION LOOKUP

function [I,pow10] = sub_find(opt,str,cstr,flag)
%[I,pow10] = sub_find(opt,str,cstr,flag_prefix)
%returns the index where str is found in cellstring array cstr
%if not found, try looking for str + prefixes

if ~ischar(str); warn('expecting a string input',1); end
if ~iscellstr(cstr); warn('expecting a cellstring input',1); end

if nargin<4; flag = true; end   %assume we want to look for prefixes
I = []; pow10 = 0;              %default values if no match found
if ~numel(str); return; end

%first pass, is a member with no prefix
tf = ismember(cstr,str);
if any(tf); I = find(tf); pow10 = 0; end

%second pass, is a member with prefix
if ~any(tf) && flag && numel(str)>1;
    tf1 = ismember(opt.prefix.symb,str(1));
    tf2 = ismember(cstr,str(2:end));
    if any(tf1) && any(tf2)
        I = find(tf2);
        pow10 = opt.prefix.powr(tf1);
    end
end

%check I doesn't have multiple hits
switch numel(I)
    case 0; %do nothing, this was not found
    case 1; %do nothing, found one match
    otherwise;
        fprintf(1,'finding "%s"\nin:\n',str);
        disp(char(cstr));
        warn('multiple matches found',1);
end
end

function ud = sub_base_unit(opt,str)
%look up base units for expression "str"
%out = sub_base_unit(opt,str)

targ = opt.unit; %target is the unit set
ud.symb = str;

%split into individual units
[u,p] = sub_parse_units(str);
n = u; s = 1; %scaling

%find each unit in the list of symbols
for j=1:numel(u)
    
    %match on symbol, with or without prefix
    [I,pow10] = sub_find(opt,u{j},targ.symb);
    %match on name, with or without prefix
    if ~numel(I); [I,pow10] = sub_find(opt,u{j},targ.name); end
    %match on name with plural "s" after. Warn if this is done
    if ~numel(I);
        [I,pow10] = sub_find(opt,u{j}(1:(end-1)),targ.name);
        warning([mfilename ':BadUnit'],...
            'looking for unit %s, presume you meant %s\n',...
            u{j},targ.name{I});
        u{j} = regexprep(u{j},'s$','');
    end
    
    %if not found, warn but continue
    if ~numel(I); warn(['can''t find unit "' u{j} '"'],0);  keyboard; end
    
    n{j} = targ.expr{I};
    s = s*(((10^pow10)*targ.size(I))^p(j));
    
end

%express as string of dimensions
ud.expr = sub_replace_units(opt,str,u,n,opt.debug); %what happens if we don't use
ud.expr = sub_simplify(ud.expr);

%this ?

ud.dims = '';
for i=1:numel(n);
    [a,b] = sub_parse_units(n{i});
    for j=1:numel(a)
        ud.dims =[ud.dims a{j} num2str(b(j)*p(i)) '.'];
    end
end
ud.size = s;
end

function [str,s] = sub_dimensions(opt,str,s,vals,flag_dims_only,flag_verb)
%[str,s] = sub_dimensions(opt,str,s,vals,flag_dims_only,flag_verb)
%converts the unit in str to basic dimensions and a scale

if nargin<5; flag_dims_only = false; end
if nargin<6; flag_verb = false; end

%do we allow reference to other "base" units, or just "kind"s of unit.
switch vals
    case 'base';    full.symb = opt.unit.symb;
        full.expr = opt.unit.expr;
        full.size = opt.unit.size;
        full.name = opt.unit.name;
    case 'kind';    full.symb = opt.kind.symb;
        full.expr = opt.kind.expr;
        full.size = opt.kind.size;
        full.name = opt.kind.name;
    otherwise
        dbstack
        disp(vals)
        error('should be ''base'' or ''kind'' here');
end
if iscell(full.size); keyboard; end

%do we seek to reduce to dimensions or just defined kinds
if flag_dims_only;  targ = opt.dims;
else                targ = opt.kind;
end

flag_done = false;
while ~flag_done;
    if flag_verb; fprintf(1,'=== starting new pass, unit is currently %s scale %f\n',str,s); end;
    flag_done = true;
    [u,p] = sub_parse_units(str);
    n = u;
    %go through each term and replace with a base unit
    %do these tests in order, ending the moment you get a hit on any test
    for j=1:numel(u)
        %if flag_verb; fprintf(1,'unit %s term %s\n',str,u{j}); end
        
        %test 1 & 2: base symbols, poss with prefix
        [I,pow10] = sub_find(opt,u{j},targ.symb);
        if numel(I);
            s = s*(targ.size(I).*10^pow10)^p(j);
            n{j} = targ.symb{I};
            continue;
        end
        
        %test 3 & 4 : full symbols, poss with prefix
        [I,pow10] = sub_find(opt,u{j},full.symb);
        if numel(I);
            s = s*(full.size(I).*10^pow10)^p(j);
            n{j} = full.expr{I};
            continue
        end
        
        %if we are here, none of the above hit
        warn(['unit ' u{j} ' not found']);
    end
    
    %if any terms changed, replace them, and redo loop above
    if ~isequal(u,n);
        flag_done = false; %redo loop above
        str = sub_replace_units(opt,str,u,n,flag_verb);
        str = sub_simplify(str,flag_verb);
        if flag_verb; fprintf(1,'unit is %s x %f at end of pass\n',str,s); end
        if s==0; keyboard; end
    end
end
str = sub_simplify(str,flag_verb);
end

function s = sub_replace_units(~,str,u,n,flag_verb)
%replaces units u with n in string str
%s = sub_replace_units(opt,str,u,n,flag_verb)

if nargin<5; flag_verb = false; end

% bugfix 20140113 compound unit powers not being handled properly

%split input into name/power
[v,p] = sub_parse_units(str);

%20150129 BUG only allow one substitution per term, do not let them cascade
%go through each sustitution, but only perform one substitution per term
for i=1:numel(v)
    old = v{i};
    for j=1:numel(u)
        new = regexprep(old,u{j},n{j});
        if ~isequal(new,old);
            v{i} = new;
            break
        end
    end
end

%create new string from (altered) v,p
s = '';
for i=1:numel(v)
    %for each term of v, split onto subunits, then apply power to each
    [a,b] = sub_parse_units(v{i});
    for j=1:numel(a)
        s = [s a{j} num2str(b(j)*p(i)) '.']; %#ok<AGROW>
    end
end
s(end) = '';

if flag_verb; fprintf(1,'replaced units "%s"->"%s"\n',str,s); end

end

function [n,pow] = sub_parse_units(str)
%splits string str into units n and powers pow
%[n,pow] = sub_parse_units(str)

%split into terms
%each term may start with / . * or whitespace \s
%then a number of letters
%then some other stuff until the next match
%TODO get non-greedy flag to work with regexp, so get all terms with one
%call.

if isempty(str); n = {}; pow = []; return; end

%if string starts with 1/unit, replace with /unit
pat = '^1/[a-z|A-Z]'; if regexp(str,pat); str(1) = ''; end

pat = '[\.|\*|/|\s+]?[a-z|A-Z]+';
beg = regexp(str,pat,'start');
if beg(1)~=1; error('unit should start with a letter or /'); end
bak = unique([beg(2:end)-1 numel(str)]);
for i=1:numel(beg);
    terms{i} = str(beg(i):bak(i)); %#ok<AGROW>
end

%if a term starts with /, invert the power of all following terms
%e.g. kg/ms is kg1m-1s-1.
pat = '^/';
pow = double(cellfun('isempty',regexp(terms,pat)));
pow(pow==0) = -1;
pow = cumprod(pow);
terms = regexprep(terms,pat,'');

%remove leading . and * from each term
pat = '^[\.|\*|\s+]';
terms = regexprep(terms,pat,'');

%parse and remove the units in each term
pat = '^[a-z|A-Z]+';
n = regexp(terms,pat,'match','once');
terms = regexprep(terms,pat,'');

%each term should now just be a power.
%remove starting ^ before powers
terms = regexprep(terms,'^\^','');

%remove ending .
terms = regexprep(terms,'\.$','');

%parse powers term-by-term to catch fractions (/)
for i=1:numel(terms)
    if isempty(terms{i}); val(i) = 1; continue; end %#ok<AGROW>
    if regexp(terms{i},'/')
        a = regexp(terms{i},'(?<num>.+)/(?<den>.+)','names');
        val(i) = str2double(a.num)/str2double(a.den); %#ok<AGROW>
    else
        val(i) = str2double(terms{i}); %#ok<AGROW>
    end
end
pow = pow .* val;

%consolidate multiples of the same unit
if numel(unique(n))<numel(n)
    [u,~,c] = unique(n);
    p = [];
    for i=1:numel(u);
        p(i) = sum(pow(c==i)); %#ok<AGROW>
    end
    n = u;
    pow = p;
end

if any(isnan(pow)); keyboard; end
end

function strud = sub_simplify(str,flag_verb)
%simplifies the units of str, removing zero-power units
%str = sub_simplify(str,flag_verb)

if nargin<2; flag_verb = false; end

%if cell array, recursively call self
if iscellstr(str)
    for i=1:numel(str)
        strud{i} = sub_simplify(str{i},flag_verb); %#ok<AGROW>
    end
    return
end

%simplify expressions
[n,pow] = sub_parse_units(str);
%remove any with zero power and consolidate into unique names
keep = pow~=0;
pow = pow(keep);
n = n(keep);
%create string out
strud = '';
if ~isempty(n)
    for u = unique(n)
        [~,loc] = ismember(u,n);
        p = sum(pow(loc));
        strud = [strud u{1} num2str(p) '.']; %#ok<AGROW>
    end
end
if flag_verb; fprintf(1,'simplify "%s" -> "%s"\n',str,strud); end

end

function opt = sub_temp_convert(opt)
%check input is purely temperature
%unit check has been performed by this point, so know input/output match

%check unit is just temperature
pat = '^T[\d|\.]+$';
if isempty(regexp(opt.in.dims,pat,'once')); error('units must be just temperatures when using abstemp option'); end
if isempty(regexp(opt.ud.dims,pat,'once')); error('units must be just temperatures when using abstemp option'); end

%check unit is a single temperature (no C/F C*K and so on)
[nin,pin] = sub_parse_units(opt.in.symb);
[nud,pud] = sub_parse_units(opt.ud.symb); %#ok<NASGU>

if ~isequal(pin,1); error('can only use the "abstemp" option with temperatures to the power of 1'); end
if numel(nin)>1; error('absolute temperature conversion can only accept a single unit at a time'); end
if numel(nud)>1; error('absolute temperature conversion can only accept a single unit at a time'); end
nin = nin{1}; nud = nud{1};

%express all in Celcius
switch lower(nin);
    case {'c','celcius','degc'};   tin = opt.amt;
    case {'k','kelvin'};    tin = opt.con.K2C(opt.amt);
    case {'f','farenheit'}; tin = opt.con.F2C(opt.amt);
    case {'ra','rankine'};  tin = opt.con.R2C(opt.amt);
    otherwise; error('unknown temperature unit');
end

%convert to output units
switch lower(nud)
    case {'c','celcius'};   tud = tin;
    case {'k','kelvin'};    tud = opt.con.C2K(tin);
    case {'f','farenheit'}; tud = opt.con.C2F(tin);
    case {'ra','rankine'};  tud = opt.con.C2R(tin);
    otherwise; error('unknown temperature unit');
end

opt.in.size = opt.amt;
opt.ud.size = tud;
end

%% OUTPUT

function sub_check_units(opt)
%check that input/output units are the same
[nin,pin] = sub_parse_units(opt.in.dims);
[nud,pud] = sub_parse_units(opt.ud.dims);
if ~isequal(nin,nud)
    disp('dimensions in:');   sub_disp_powers(opt,nin,pin); disp(' ');
    disp('dimensions out:');  sub_disp_powers(opt,nud,pud);
    error('dimensions in/out do not match');
end

if isempty(nin); return; end


%warn or error on power mismatch
err = abs(pin - pud);
worst = max(err);

tol = [1E-6 1E-3];

%if all within 1E-6, don't care
if worst<tol(1); return; end

%if any worse than 1E-3, error
if worst>tol(2);
    disp('dimensions in'); sub_disp_powers(opt,nin,pin); disp(' ');
    disp('dimensions out');sub_disp_powers(opt,nud,pud); disp(' ');
    error(['unit powers do not match to within ' num2str(tol(2)) ', worst error ' num2str(worst)]);
end

%warn about slight mismatch, but continue
if worst>=tol(1) && worst<=tol(2)
    warning(['unit powers do not match, worst error ' num2str(worst)])
end

    function sub_disp_powers(opt,n,p)
        for i=1:numel(n);
            [~,loc] = ismember(n{i},opt.kind.symb);
            fprintf('%s^%g\n',opt.kind.name{loc},p(i));
        end
    end
end

function str = sub_pretty_number(opt,num)

%engineering format, if requested
%to nearest 1E3 below
if opt.eng
    pow = floor(log10(num)/3)*3;
    str = sprintf('%gE%d',num/(10^pow),pow);
    return
end

%is it a round fraction ?
[n,d] = rat(num);
if d==1 && n>0; str = sprintf('%d',n);  return; end
if n==1;        str = sprintf('%d/%d',n,d); return; end

%is it a pi fraction ?
[nr,dr] = rat(num/pi);
if dr==1 && nr==1;      str = 'pi';                     return; end
if dr==1 && nr>1;       str = sprintf('%gpi',nr);       return; end
if nr==1;               str = sprintf('pi/%d',dr);      return; end

%is it a pi multiple ?
[nr,dr] = rat(num*pi);
if nr==1 && dr==1; str = '1/pi';         return; end
if dr==1 && nr>1 ; str = sprintf('%d/pi',nr);   return; end
if nr==1;   str = sprintf('1/%dpi',dr); return; end

%default, g format (if none of the others applied
str = sprintf('%G',num);

end

function sub_show_units(opt)
if opt.abstemp;
    for i=1:numel(opt.in.size);
        fprintf(1,'%s %s\t=\t%s %s\n',...
            sub_pretty_number(opt,opt.in.size(i)),...
            opt.in.symb,...
            sub_pretty_number(opt,opt.ud.size(i)),...
            opt.ud.symb);
    end
    return
end

v = opt.conv.*opt.amt;
for i=1:numel(opt.amt);
    fprintf(1,'%s %s\t=\t%s %s\n',...
        sub_pretty_number(opt,opt.amt(i)),...
        opt.in.symb,...
        sub_pretty_number(opt,v(i)),...
        opt.ud.symb);
end
disp('-');

fprintf(1,'%s\t=\t%s\tx %s\n',opt.in.symb,opt.ud.symb,sub_pretty_number(opt,opt.conv));
fprintf(1,'%s\t=\t%s\tx %s\n',opt.ud.symb,opt.in.symb,sub_pretty_number(opt,1/opt.conv));
disp('-');
sub_show_dimensions(opt,opt.in,'from');
disp('-');
sub_show_dimensions(opt,opt.ud,'to');

end

function sub_show_dimensions(opt,in,str)

if nargin<3; str = ' '; end
fprintf(1,'%-6s: %s\n',str,in.symb);
fprintf(1,'scale : x%s of base\n',sub_pretty_number(opt,in.size));

[n,p] = sub_parse_units(in.dims);
ns = opt.kind.name(opt.kind.isdim);
us = opt.kind.unit(opt.kind.isdim);
ss = opt.kind.symb(opt.kind.isdim);
ps = zeros(size(ns));
for ii=1:numel(ns);
    [tf,loc] = ismember(ss{ii},n);
    if tf; ps(ii) = p(loc); end
end

%print numerator
if ~opt.frac;
    str1 = sub_np2str(opt,ns,ps);
    str2 = sub_np2str(opt,us,ps);
    fprintf(1,'%s : %s\n','base ',str2);
    fprintf(1,'%s : %s\n','kind ',str1);
    
else
    strnum = sub_np2str(opt,us(ps>0),ps(ps>0));
    strden = sub_np2str(opt,us(ps<0),-ps(ps<0));
    fprintf('%-6s: %s /\n','base',strnum);
    %len = max([numel(strnum) numel(strden)]);
    %fprintf('%6s: %s\n','units',repmat('-',[1 len]));
    fprintf('%-6s  %s\n','     ',strden);
    
    strnum = sub_np2str(opt,ns(ps>0),ps(ps>0));
    strden = sub_np2str(opt,ns(ps<0),-ps(ps<0));
    fprintf('%-6s: %s /\n','kind',strnum);
    %len = max([numel(strnum) numel(strden)]);
    %fprintf('%6s: %s-\n','kind',repmat('-',[1 len]));
    fprintf('%-6s  %s\n','    ',strden);
    
end


function str = sub_np2str(opt,n,p)
str = '';
for i=1:numel(n);
    if ~p(i); continue; end
    if p(i)==1;
        str = [str sprintf(' x %s',n{i})]; %#ok<AGROW>
    else
        str = [str sprintf(' x %s^%s',n{i},sub_pretty_number(opt,p(i)))]; %#ok<AGROW>
    end
end
if numel(str); str = str(4:end); end
end
end

function warn(str,halt)
%this is a catcher function for warning to allow debugging
%in everyday use, should drop an error and exit
if nargin<2; halt = 1; end
if isempty(str); str = ' '; end

switch halt;
    case -1; error(str);
    case 0;
        warning([mfilename ':Problem'],str);
    otherwise
        warning([mfilename ':Problem'],str);
        disp('type "dbquit" to exit from this warning');
        keyboard;
end
end

function sub_report(opt)
disp('==================')
disp('== UNIT REPORT ===');
disp('==================');
fprintf('These are the inputs recognised by the function %s\n',mfilename);

disp('================');
disp('=== OPTIONS ====');
disp('================');
disp(char(setxor(fieldnames(opt),{'amt','prefix','con','kind','dims','unit','in','ud','newunit'})));
disp(' ');

disp('================')
disp('=== PREFIXES ===');
disp('================')
for i=1:numel(opt.prefix.symb)
    fprintf('%c 10E%+d\n',...
        opt.prefix.symb(i),...
        opt.prefix.powr(i));
end
disp('Each term in a unit expression is allowed a single prefix')
disp('For example W, kW, MW, GW are all multiples of watts')
disp(' ');

disp('======================');
disp('=== KINDS OF UNITS ===');
disp('======================');
fprintf('%12s %5s %5s %6s %5s %10s\n',...
    'name','symb','expr','unit','dim?','dims');
for i=1:numel(opt.kind.name)
    fprintf('%12s %5s %5s %6s %5d %10s\n',...
        opt.kind.name{i},...
        opt.kind.symb{i},...
        opt.kind.expr{i},...
        opt.kind.unit{i},...
        opt.kind.isdim(i),...
        opt.kind.dims{i});
end
disp('dim? shows if this is a basic dimension')
disp('for example "length" is a basic dimension, but "area" is not, as area has dimension length^2')
disp('all units in and out are reduced to basic dimensions to check if dimensions match');
disp(' ');
disp('========================')
disp('=== UNITS (BY NAME) ====');
disp('========================')
fprintf('%16s %5s %15s\n',...
    'NAME','SYMB','EXPRESSION');
%sort units by name
[~,order] = sort(lower(opt.unit.name));

for i=1:numel(opt.unit.name)
    fprintf('%16s %5s %8s x %s\n',...
        opt.unit.name{order(i)},...
        opt.unit.symb{order(i)},...
        sub_pretty_number(opt,opt.unit.size(order(i))),...
        opt.unit.expr{order(i)});
end

disp(' ');
disp('=========================')
disp('=== UNITS (BY SYMBOL) ===')
disp('=========================')
fprintf('%16s %5s %15s\n',...
    'NAME','SYMB','EXPRESSION');
%sort units by name
[~,order] = sort(lower(opt.unit.symb));
for i=1:numel(opt.unit.name)
    fprintf('%16s %5s %8s x %s\n',...
        opt.unit.name{order(i)},...
        opt.unit.symb{order(i)},...
        sub_pretty_number(opt,opt.unit.size(order(i))),...
        opt.unit.expr{order(i)});
end

disp(' ');
disp('your input/output units must all be from the "SYMB" or "NAME" column above');
disp('For example ''V'' will be interpreted as Volts, not Volume - because volts is a unit, and volume is a kind of unit')
disp('A volume unit would be m3 or cc')
disp(' ');
disp('terms must be separated by one of *-./ e.g. N-m or m/s')
disp('each term may have zero or one prefixes, e.g. kN/m or kg/cm3')
disp('each term may have a zero or one power, a number, e.g. kg/m3 or lbf*in-2.')
disp('powers may be fractions, e.g. kg^1/2 or mm1/4');
disp(' ');
disp('The expressions in this table are all references to the "symb" column of the kinds table')
end

function sub_list(opt)
disp('=========================');
disp('=== LIST OF ALL UNITS ===');
disp('=========================');
a = [opt.unit.name opt.unit.symb];
a = unique(a);
disp(char(a));
disp('=========================');
end

%% DEVLOG
%2014.07.01 added "sec" as unit relating back to "s"
%2014.07.08 added eng,exact,rat options
%2014.07.08 added display of pi and 1/pi multiples
%2014.07.09 bugfix if third input is mutiple values, only show first
%2014.07.22 accept units in any order
%TODO extra flag for how to show units
%TODO add vrms,irms as functions of power
%2014.08.13 removed exact,rat options
%2014.08.13 added 'abstemp' option
%2014.08.13 fixed bug where preample+symbolname was not parsing (e.g. kohm)
%2015.01.26 DONE add units on the fly (as structures)
%2015.01.26 DONE check new units correct
%2015.01.26 DONE fixed "rankine" changed symbol "R"->"Ra" to resolve clash with "ohm"
%2015.01.26 DROP infosize change symbol s->i
%2015.01.25 DONE defining rpm as rev/min goes wrong, but rev/s okay
%2015.01.28 DONE allow for whitespace as separator between units
%2014.01.28 DONE "frac" option to show output dimesions as fractions
%2015.01.28 DONE catch known unit with "s" on the end

%% DEVNOTES
%There are two tables: "kind" and "unit"
%kind is classes (e.g. mass)
%unit is different units within each class (e.g. pounds, grammes, etc)
%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% DIMENSION CHECKS %%%
%%%%%%%%%%%%%%%%%%%%%%%%
%A subset of the kinds are dimensions - irreducible units.
%when doing the dimension check, all units are reduced to their basic
%dimensions. this is done by iteratively replacing each unit with the expression for
%that unit (and appropriate scaling) until only dimensions remain.
%obviously, the unit expression have to reduce to basic dimensions, or this
%will fail.
%example:
%the kind of unit "force" is not a dimension.
%it is replaced by the expression for force m.a (mass * acceleration)
%this is repeated:
%force ->
%mass * acceleration ->
%mass * velocity / time ->
%mass * length / time^2
%these are all base dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADDING NEW UNITS %%%
%%%%%%%%%%%%%%%%%%%%%%%%
%new units go in sub_base_data
%new kinds go in sub_data_data
%
%All unit symbols should be unique.
%Also, if a unit matches a prefix/unit pair, the unit will be used:
%e.g. if you define a new unit "Hectare" with symbol "HA"
%all references to hecto-amp will be lost, because it also has the symbol "HA"
%
%avoid introducing circular references in the expressions for a unit.