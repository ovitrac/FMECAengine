function [gaussianpeak,outsum,outsumbaseline] = monotonepeakfit(p,varargin)
%MONOTONEPEAKFIT fits [x y] as a sum of Gaussian/Lorentzian peaks with positions previsously identified with monotonepeak
% Note that two methods of fit are always applied with prescribed centers (method 1) and movable centers (method 2)
% Results of both methods are returned.
% The simplified Gaussian kernel used in this function verifies: 
%                 area: 0.6006*sqrt(pi)*width
%   standard deviation: 0.6006/sqrt(2)*width
% The corresponding Lorenzian kernel shares the same value at half height
%
% SYNTAX: gaussianpeak = monotonepeakfit(p,'x',x,'y',y,'property1',value1...,'keyword')
%         [gaussianpeak,outsum,outsumwithbaseline] = monotonepeakfit(...)
%
% LIST OF PAIR PROPERTY/VALUE
%                p: npx1 structure array created with MONOTONEPEAK (with np = peak number)
%                x: vector of x values (default = [])
%                y: mx1 vector of y values
%      significant: percent to keep signifiant peaks (default = .95)  % search the first peaks >95% of weight
%     minpointsinbaseline: (default 5)
%         shiftmax: maximum shift tolerated with method 2 (default = []), scalar or 1xnp vector to customize shifts for each peak
%      shiftbuffer: shift tolerance buffer (default = []), scalar quantity or or 1xnp vector to customize shifts for each peak
%shiftpenaltyscale: shift Langrange multiplier (default = []), scalar or 1xnp vector to enforce constraints differently for each peak
%         widthmax: maximum width (default = Inf), scalar or 1xnp vector to customize shifts for each peak
%      widthbuffer: shift tolerance buffer (default = []), scalar or 1xnp vector to customize shifts for each peak
%widthpenaltyscale: width Langrange multiplier (default = []), scalar or 1xnp vector to enforce constraints differently for each peak
%          preject: rejection probability (percentile) for baseline detection (default = 5);
%
% KEYWORDS
%       'baseline': fit and remove linear baseline 
%           'sort': to sort in an descend order of weigth attributed to gaussian peaks
%     'lorentzian': to fit lineshape with lorentzian peak
%      'endforced': to add last point of x to identify baseline
%      'keeporder': force Gaussians/Lorentzian to have increasing positions
%  keepfittedorder: keeps the peak order as fitted
% keepinitialorder: keeps the initial peak order
%     'independent': start strategy 2 with the same guess as strategy 1 instead of using the width of strategy 1
% 'dispconstraints': print when contraints are introduced
%
% OUTPUTS (with the following notations: np = number of peaks, 2 = 2 fitting models)
%
%    gaussianpeak = np x 2 structure array with fields
%              rank: rank of gaussian peak (descendant order by weight)
%            weight: weight of gaussian peak 
%    relativeweight:relative weight of gaussian peak compaired to sum of weight
%             width: width of gaussian peak
%          position: position or center of gaussian peak
%          baseline: polyval function to define baseline of gaussian peak
%        isbaseline: index/point to define baseline of gaussian peak
%         xbaseline: point based on x to define baseline of gaussian peak
%            kernel: kernel function of gaussian peak
%            nonneg: replace G\y by lsqnonneg(G,y) for the amplitude of peaks if true (default = false)
%                r2: coefficient of determination of gaussian peak
%             r2max: 2 maxima of coefficient of determination of gaussian peak
%             sigma: standard deviation of gaussian peak
%        expansion1: kernel function of all gaussian peaks without baseline (x, n = peak number, k = fitting models)
%        expansion2: kernel function of all gaussian peaks with baseline (x, n = peak number, k = fitting models)
%            window: is a structure with fields
%                   center: center of window 
%                    width: width of window (x(end)-x(1))
%                widthapha: width of window calculated with alpha
%                    alpha: probability of rejection in window (default=0.05)
%        nsignificantpeaks: number of significant peaks
%         significantproba: percentage to  keep signifiant peaks (variable 'significant', default = .95)
% 
%   outsum is an anonymous function (x,norder,imodel) giving the serial expansion
%               Sum i =1..norder weight(i) * kernel_i(x)      for imodel = 1 or 2
%   outsumwithbaseline is an anonymous function (x,norder,imodel) giving the serial expansion
%               (Sum i =1..norder weight(i) * kernel_i(x)) + baseline(x)       for imodel = 1 or 2
%
% See also monotonepeak, monotone
%
% BASIC example
%   g = @(x,position,width)  exp(-((x-position)./(0.6006.*width)) .^2)
%   x = linspace(0,1000,1000)';
%   y = 0.5*g(x,700,30) + 2*g(x,830,100);
%   p = monotonepeakfit; p.center = [690 850]; p.width = [10 20];
%   q = monotonepeakfit(p);
%   [q,~,model] = monotonepeakfit(p,'x',x,'y',y);
%   hp = plot(x,[y model(x,2,1) model(x,2,2)]); legend(hp,{'original' 'fit 1' 'fit with free center'},2)
%   
%
% ADVANCED EXAMPLE
% [dbpur,dbxpur] = nmrloadbase; 
% mol = {'Erucamide'}; % substance test Erucamide
% [~,isubs] = intersect(dbxpur.commonname,mol); % find index of Erucamide in databse
% ROI = dbpur.(mol{1}).gates; % extract RoI (nx4 array: ppmmin, ppmmax, buffer and weight)
% valid = ((dbxpur.ppm>=ROI(1,1))&(dbxpur.ppm<=ROI(1,2))); % take the 1st considered ROI
% x = dbxpur.ppm(valid);          
% y = dbxpur.I0(valid,isubs);                
% p = monotonepeak('x',x,'y',y,'array'); % find peaks
% [gaussianpeak,~,model] = monotonepeakfit(p,'x',x,'y',y,'baseline','sort'); % fitting
% subplot(211)
% hp = plot(x,[y model(x,size(gaussianpeak,1),1) model(x,size(gaussianpeak,1),2)]);
% hl = legend('measured','deconvoluted (s)','deconvoluted (p,s)');
% subplot(212)
% hp = plot(x,[y model(x,3,1) model(x,3,2)]);
%
%
% INRA\MS 2.1 - 24/03/2013 - Olivier Vitrac - rev. 07/07/2018
%
%
% TODO LIST
%   16/05/13 help on additional parameters, fix r2
%
% Revision history
% 25/03/13 add help and example
% 26/03/13 major update, remove cells, add outsum, add r2all, fix r2
% 27/03/13 add lorentzian peaks and relativeweight
%          fix help
% 29/03/13 add gaussianpeak.window for exact description of multiplet (center, width...)
% 09/04/13 add variable 'significant' to define rejection percent 
%          modify outsum = model without baseline
%                 outsumbaseline = model with baseline 
% 11/04/13 rename variables and add window.nsignificantpeaks,window.significantproba
% 16/04/13 add keyword 'endforced' to add last point of x to identify baseline
% 14/05/13 fix outputs, major help update
% 15/05/13 update help
% 23/11/13 accept mfilt=0
% 30/10/15 implements shiftmax as vector, add shiftbuffer and shiftpenaltyscale
% 31/10/15 fix for x in decreasing order (H(x) has been modifidied consequently so that dx>0 in all cases)
% 31/10/15 add widthmax, widthbuffer, widthpenaltyscale, independent, keepinitialorder
% 01/11/15 harmonization of constraints handling
% 10/11/15 add dispconstraints
% 11/11/15 change argcheck behavior to expand structures while protecting options ('nostructaxpand')
% 14/11/15 add area field
% 06/12/15 replace keepinitialorder by keepfittedorder, implements a true keepinitialorder based on nearestpointunique
% 29/06/18 add nonneg -> rank check
% 04/07/18 add widthmin
% 07/07/18 fix when all cumweights are NaN

%% default
keyword = {'baseline','sort','lorentzian','endforced','keeporder','keepinitialorder','keepfittedorder','independent','dispconstraints'};
options = struct('display','iter','FunValCheck','on','MaxIter',1e3,'TolFun',1e-6,'TolX',1e-6);
default = struct('x',[],'y',[],'significant',.95,'minpointsinbaseline',5,...
            'shiftmax',[],'shiftbuffer',[],'shiftpenaltyscale',[],...
            'widthmax',Inf,'widthmin',-Inf,'widthbuffer',[],'widthpenaltyscale',[],...
            'preject',5,'nonneg',false);
minumunrequiredfields = {'center','width','ibase'};
requiredfields = {'tail' 'wall' 'height' 'center' 'istart' 'start' 'stop' 'istop' 'width','ratioheight','ratiowidth','ibase'};
maxpeaksforconstraints = 7;

%% Kernels
% simplified Gaussian kernel
%   area: 0.6006*sqrt(pi)*width
%   standard deviation: 0.6006/sqrt(2)*width
% Lorenzian kernel shares the same value at half height
gaussiankernel = @(x,position,width)  exp(-((x-position)./(0.6006.*width)) .^2);
lorentziankernel =  @(x,position,width) 1./(1+((x-position)./(0.5.*width)).^2);

% argcheck
if nargin<1, gaussianpeak = cell2struct(repmat({[]},1,length(minumunrequiredfields)),minumunrequiredfields,2); return, end
if isempty(p), p = monotonepeak(varargin{:}); end
if ~isstruct(p) && all(cellfun(@(f) isfield(p,f),minumunrequiredfields)), error('p must be created with monotonepeak'), end
ismonotonepeakobject = all(cellfun(@(f) isfield(p,f),requiredfields));
if ~ismonotonepeakobject, dispf('WARNING:: object p has not been created with monotonepeak'), end
if length(p)==1 && length(p.center)>1, p = struct2structtab(p); end
n = length(p);
% retrieve options with 'nostructexpand' and bring the reminaing to normal argcheck (implemented on 11/11/15)
[otmp,remain] = argcheck(varargin,struct('options',options),keyword,'nostructexpand');
o = argcheck(remain,default);
o = argcheck(otmp,o,'','keep');
m = length(o.y);
if m == 0, error('y is empty'), end
if numel(o.y)~=m, error('y must be a vector'), end
if isempty(o.x), o.x = (1:m)'; end
if numel(o.x)~=m, error('x and y must be of the same size'); end
o.x = o.x(:); o.y = o.y(:);
sce = var(o.y)*(m-1);
penaltyscale = sqrt(sce)/n;
xrange = [min(o.x) max(o.x)]; xwidth = diff(xrange);
dx = abs(median(diff(o.x))/4);
H = @(x) 1/2 * ( 1 + tanh(x/dx) ); % heaviside
% keyword 'lorentzian' if fitting with lorentzian peaks
if o.lorentzian, gaussiankernel = lorentziankernel; end

% Constraints control (shift and width), note that they are used only for the second fitting strategy
if isempty(o.shiftmax), o.shiftmax = min([p.width]); end, o.shiftmax = abs(o.shiftmax);
if isempty(o.shiftbuffer), o.shiftbuffer = dx; end
if isempty(o.shiftpenaltyscale), o.shiftpenaltyscale = penaltyscale; end
if isempty(o.widthmax), o.widthmax = o.shiftmax; end
if isempty(o.widthmin), o.widthmin = -Inf(size(o.widthmax)); end
if isempty(o.widthbuffer), o.widthbuffer = o.shiftbuffer; end
if isempty(o.widthpenaltyscale), o.widthpenaltyscale = o.shiftpenaltyscale; end
for prop = {'shiftmax' 'shiftbuffer' 'shiftpenaltyscale' 'widthmax' 'widthmin' 'widthbuffer' 'widthpenaltyscale'}
    o.(prop{1}) = argpad(o.(prop{1}),n);
    o.(prop{1}) = o.(prop{1})(:)';
end
xbuffer = o.shiftbuffer; Hshift = @(x) 1/2 * ( 1 + tanh(x./xbuffer) ); % heaviside for shift constraint
xbuffer = o.widthbuffer; Hwidth = @(x) 1/2 * ( 1 + tanh(x./xbuffer) ); % heaviside for width constraint


% remove baseline
if o.baseline
    isbaseline = true(m,1);
    for i=1:n
        isbaseline( (o.x>=p(i).start) & (o.x<=p(i).stop) ) = false;
    end
    if length(find(isbaseline))<o.minpointsinbaseline
        isbaseline([1 end])=true;
        isbaseline([p.ibase]) = true;
    end
    if o.endforced
        fluctuations = o.y;
        fluctuationsreject = prctile(fluctuations,o.preject);
        isbaseline(1:find(fluctuations>fluctuationsreject,1,'first')) = true;
        isbaseline(find(fluctuations>fluctuationsreject,1,'last'):m) = true;
    end
    b = polyfit(o.x(isbaseline),o.y(isbaseline),1);
    o.y = o.y - polyval(b,o.x);
else
    b = [0 0];
end

% Fitting
G = zeros(m,n);       % kernels (used by nested functions)
widthguess = cat(2,p.width)/4; %first guess of s
positionguess = cat(2,p.center); % first guess of p
[width,position,sigma,area]  = deal(NaN(2,n)); % 2 solutions - row-wise
[weight,cumweight,order,rang,relativeweight] = deal(NaN(n,2)); % 2 solutions - column-wise
[windowcenter,nsignificantpeaks] = deal(zeros(2,1));
critfit = zeros(2,1); % 2 solutions - column-wise
yfit    = NaN(m,2);   % 2 solutions - column-wise
% STRATEGY 1:: fit based on s only (peak positions set by [p.center])
warning off %#ok<WNOFF>
[width(1,:),critfit(1)] = fminsearch(@fitgaussianfixedpos,widthguess,o.options);
position(1,:) = positionguess;
% STRATEGY 2:: fit with positions and s free
if o.independent
    [tmp,critfit(2)] = fminsearch(@fitgaussian,[positionguess;widthguess],o.options);
else
    [tmp,critfit(2)] = fminsearch(@fitgaussian,[positionguess;width(1,:)],o.options);
end
position(2,:) = tmp(1,:); width(2,:) = tmp(2,:);

% Strategy 2 may change the positions of peaks, this section try to find the combination, which matches user expectations
% NB: the sum of Gaussian is commutative, the math itself cannot resolve the issue (constraints may help)
% try to find the combination which minimizes the number of constraint violations
if n<=maxpeaksforconstraints
    perm = fliplr(perms(1:n)); nperm = length(perm); violations = zeros(nperm,7);
    [~,p0order] = sort(positionguess); [~,porder] = sort(position(2,:)); % sorted peak positions
    widthmax = max(width(2,:));
    for i = 1:nperm
        violations(i,1) = sum(abs(porder(perm(i,:))-p0order)); % peak order
        violations(i,2) = sum((abs(position(2,perm(i,:))-positionguess))); % position deviation (to be normalized)
        violations(i,3) = sum(abs(width(1,:)-width(2,perm(i,:)))); % width deviation (to be normalized)
        violations(i,4) = sum((abs(position(2,perm(i,:))-positionguess))>o.shiftmax)/n; % to large shifts
        violations(i,5) = sum(width(2,perm(i,:))>o.widthmax)/n; % unverified peak widths
        violations(i,6) = sum(width(2,perm(i,:))>3*4*widthguess)/n; % unverified peak widths
        violations(i,7) = sum(width(2,perm(i,:))<widthguess)/n; % unverified peak widths
    end
    maxviolations = max(violations,[],1); violations(:,1:3) = violations(:,1:3)./maxviolations(ones(nperm,1),1:3);
    if widthmax>xwidth/2,  % widths are so large that positions are unusable
                              violationweights = [0  ; 0 ;1  ;0  ;1  ;1 ; 1 ]; % current weight for each rule
    elseif widthmax>xwidth/4, violationweights = [.5 ;.5 ;1 ;.5  ;1  ;1 ; 1]; % current weight for each rule
    else                      violationweights = [1.5;1  ;1  ;1  ;.3 ;.3;.3]; % current weight for each rule
    end
    [~,ibest] = min(violations * violationweights);
    position(2,:) = position(2,perm(ibest,:));
    width(2,:) = width(2,perm(ibest,:));
end

% distance
warning on %#ok<WNON>
numrank = zeros(1,2);
for k=1:2
    [yfit(:,k),weight(:,k),numrank(k)] = gaussian(o.x,position(k,:),width(k,:));
    [~,order(:,k)] = sort(weight(:,k),'descend');
    rang(order(:,k),k) = 1:n;
end

% control plot (to be removed)
% figure, plot(o.x,[o.y yfit]); hold on, plot(o.x(isbaseline),o.y(isbaseline),'r.')

% list of peaks
gaussianpeak = repmat(struct('rank',[],'weight',[],'relativeweight',[],...
                             'width',[],'sigma',[],'area',[],'position',[],'baseline',[],'xbaseline',[],'isbaseline',[],...
                             'kernel',[],'r2',NaN,'r2max',[],'numrank',numrank,...
                             'window',struct('center',[],'width',[],'widthalpha',[],'alpha',[],'nsignificantpeaks',[],'significantproba',[])),[n,2]);
alpha = 0.05; % probability of rejection
for k = 1:2
    relativeweight(:,k) = weight(:,k)/sum(weight(:,k),1);
    windowcenter(k) = position(k,:)*weight(:,k)/sum(weight(:,k),1);
    windowwidth = o.x(end)-o.x(1);
    windowwidthalpha = norminv([alpha/2 1-alpha/2],windowcenter(k),(0.6006*windowwidth)/sqrt(2));
    sigma(k,:) = (0.6006*width(k,:))/sqrt(2); % standard deviation: 0.6006/sqrt(2)*width
    area(k,:) = 0.6006*sqrt(pi)*width(k,:);   % area: 0.6006*sqrt(pi)*width
    for i=1:n
        gaussianpeak(i,k).kernel =  @(x) gaussiankernel(x,position(k,i),width(k,i));
        gaussianpeak(i,k).weight = weight(i,k);  
        gaussianpeak(i,k).relativeweight = relativeweight(i,k);
        gaussianpeak(i,k).rank = rang(i,k);
        gaussianpeak(i,k).r2max = 1 - critfit(k).^2/sce;
        if o.baseline
            gaussianpeak(i,k).baseline = @(x) polyval(b,x);
            gaussianpeak(i,k).isbaseline = isbaseline;
            gaussianpeak(i,k).xbaseline = o.x(isbaseline);
        end
        gaussianpeak(i,k).width = width(k,i);    
        gaussianpeak(i,k).position = position(k,i);        
        gaussianpeak(i,k).sigma = sigma(k,i);
        gaussianpeak(i,k).area = area(k,i);
        gaussianpeak(i,k).window.center = windowcenter(k);
        gaussianpeak(i,k).window.width = windowwidth;
        gaussianpeak(i,k).window.alpha = alpha;
        gaussianpeak(i,k).window.widthalpha = diff(windowwidthalpha);
    end
end

% sorting peaks if required
if o.sort || (nargout>1),
    
    % peak initial order
    initialorder = zeros(2,1);
    [~,ip1] = sort([gaussianpeak(:,1).position]); initialorder(1,ip1)=1:n;
    [~,ip2] = sort([gaussianpeak(:,2).position]); initialorder(2,ip2)=1:n;
    
    % sort all models
    for k=1:2
        gaussianpeak(:,k) = gaussianpeak(order(:,k),k);
        position(k,:) = position(k,order(:,k));
        width(k,:) = width(k,order(:,k));
        weight(:,k) = weight(order(:,k),k);
        sigma(k,:) = sigma(k,order(:,k));
        cumweight(:,k) = cumsum(relativeweight(order(:,k),k),1);
        if ~all(isnan(cumweight(:,k)))
            nsignificantpeaks(k) = find(cumweight(:,k)>=o.significant,1,'first');
        end
        windowcenter(k) = position(k,1:nsignificantpeaks(k))*weight(1:nsignificantpeaks(k),k)/sum(weight(1:nsignificantpeaks,k),1);       
    end
    % user override
    if o.keepinitialorder
        for k=1:2
            order = nearestpointunique(positionguess,[gaussianpeak(:,k).position]);
            gaussianpeak(:,k) = gaussianpeak(order,k);
            position(k,:) = position(k,order);
            width(k,:) = width(k,order);
            weight(:,k) = weight(order,k);
        end
    elseif o.keepfittedorder % previously o.keepinitialorder (fixed on Jan 6, 2016)
        for k=1:2
            [~,order] = sort([gaussianpeak(:,k).position]);
            gaussianpeak(:,k) = gaussianpeak(order(initialorder(k,:)),k);
            position(k,:) = position(k,order(initialorder(k,:)));
            width(k,:) = width(k,order(initialorder(k,:)));
            weight(:,k) = weight(order(initialorder(k,:)),k);
        end
     elseif o.keeporder
        [~,order] = sort([gaussianpeak(:,1).position]); gaussianpeak = gaussianpeak(order,:);
        position = position(:,order); width = width(:,order); weight = weight(order,:);        
    end
    % general expansion expression (for both models)
    expansion1 = @(x,n,k) gaussiankernel(repmat(x(:),1,n),position(k*ones(1,numel(x)),1:n),width(k*ones(1,numel(x)),1:n))*weight(1:n,k);
    expansion2 = @(x,n,k) gaussiankernel(repmat(x(:),1,n),position(k*ones(1,numel(x)),1:n),width(k*ones(1,numel(x)),1:n))*weight(1:n,k)+polyval(b,x(:));
    % use the expansion to refresh r2 and relativeweight with expansion order n
     for k=1:2
        for i=1:n
            gaussianpeak(i,k).r2 = 1 - norm(o.y-expansion2(o.x,i,k)).^2/sce;
            gaussianpeak(i,k).relativeweight =  cumweight(i,k);
            gaussianpeak(i,k).window.center = windowcenter(k);
            windowwidthalpha = norminv([alpha/2 1-alpha/2],windowcenter(k),(0.6006*windowwidth)/sqrt(2));
            gaussianpeak(i,k).window.widthalpha = diff(windowwidthalpha);
            gaussianpeak(i,k).window.significantproba = o.significant;
            gaussianpeak(i,k).window.nsignificantpeaks = nsignificantpeaks(k);
        end
     end
end

% additional outputs
if nargout>1, outsum = expansion1; end
if nargout>2, outsumbaseline = expansion2; end

%% Nested functions
    function err = fitgaussian(ps)
        WC = sum( (Hwidth( ps(2,:) - o.widthmax) + Hwidth( o.widthmin - ps(2,:)))  .*o.widthpenaltyscale );
        SC = sum((Hshift( abs(ps(1,:)-[p.center]) - o.shiftmax) ).*o.shiftpenaltyscale);
        err = norm(gaussian(o.x,ps(1,:),ps(2,:))-o.y) ...
              + norm((1-H(ps(2,:)))*penaltyscale) ... force positive variance, constraint the displacement of peaks
              + WC ... % width constraints
              + SC; % shift constraints
        if o.dispconstraints 
            if SC>0, dispf('\t -->shift constraints introduced'), end
            if WC>0, dispf('\t -->width constraints introduced'), end
        end
    end

    function err = fitgaussianfixedpos(s)
        err = norm(gaussian(o.x,[p.center],s)-o.y) + norm((1-H(s))*penaltyscale);
    end

    function [y,weights,rk] = gaussian(x,c,s)
        % base functions
        for j = 1:n
            G(:,j) = gaussiankernel(x,c(j),s(j));
        end
        % weights
        if o.nonneg && rank(G)<n
            W = lsqnonneg(G,double(o.y));
        else
            W = abs(G\o.y);
        end 
        % fit
        y = G*W;
        % additional out if required
        if nargout>1, weights = W; end
        if nargout>2, rk = rank(G); end
    end

end