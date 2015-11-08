function [RH,RHfunc] = salt2RH(salt,T)
% salt2RH gives the RH above a salt solution at temperature T (ref: http://nvlpubs.nist.gov/nistpubs/jres/81A/jresv81An1p89_A1b.pdf)
% SYNTAX:   RH = salt2RH(salt,T)
%           [RH,RHfunc] = salt2RH(...)

% Migration 2.1 - 26/10/2015 - INRA\Olivier Vitrac - rev. 27/10/2015

% Revision history
% 27/10/2015 list all salts, implements recursion on both salt and T

% Definitions
persistent dbsalt
source = struct(...
    'path',find_path_toolbox('migration'),...
    'file','RH_above_salt_solutions.ods',...
    'sheetname','data');

% load data
if ~exist('dbsalt','var') || isempty(dbsalt)
    dbsalt = loadodsprefetch(fullfile(source.path,source.file),source);
end

% arg check
if nargin<1, dispf('\nSALT22RH:: list of salts'), dispf('\t%s\n',dbsalt.salt{:}), return, end
if nargin<2, T = 25; end
if iscellstr(salt)
    salt = strtrim(salt);
    nsalt = numel(salt);
    if nsalt == numel(T)
        RH = NaN(size(salt));
        for i=1:nsalt,
            if ~ismember(lower(salt{i}),{'none' 'dummy'}) && ~isnan(T(i))
                RH(i)=salt2RH(salt{i},T(i));
            end
        end
    else
        RH = cellfun(@(s) salt2RH(s,T),salt);
    end
    return
end

% search data
i = find(~cellfun(@isempty,regexpi(dbsalt.salt,salt))); ni = length(i);
if ni>1 
    warning('\nSALT2RH:: %d matches found',ni), dispf('\t%s\n',dbsalt.salt{i})
    if nargout>1, RH = []; end
    if nargout>1, RHfunc = @() []; end
    return
elseif ni==0
    error('no match for ''%s''',salt)
end

% calculation
currentsalt = dbsalt.salt{i};
A = [dbsalt.A3(i) dbsalt.A2(i) dbsalt.A1(i) dbsalt.A0(i)];
Tmin = dbsalt.TminC(i);
Tmax = dbsalt.TmaxC(i);
RH = polval(T);
if nargout>1, RHfunc = @polval; end

% Nested function
    function RH = polval(T)
        if (T<Tmin) || (T>Tmax)
            dispf('SALT2RH::%s: current temperature (%0.5goC) is out bounds [%0.4g, %0.4g]oC',currentsalt,T,Tmin,Tmax)
        end
        RH = polyval(A,T);
    end
end