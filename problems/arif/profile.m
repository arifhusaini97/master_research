function profileSpec = profile()
%PROFILE Sweep definitions and evaluators for Arif.

modelFns = model();
vals = values();

%No 2nd solution for M=0
mValues = vals.sweep.mValues;
lambdaVals = applyLambdaExclusions(vals.sweep.lambdaVals, vals.sweep.excludeLambdaRanges);

profileSpec = struct();
profileSpec.createBaseHandles = @(cfg) buildHandles(cfg, modelFns);
profileSpec.metricEvaluators  = @(cfg) buildMetricEvaluators(cfg, modelFns);
profileSpec.figureConfigs     = @(cfg, metrics) buildFigureConfigs(cfg, metrics);
profileSpec.sweeps            = @(cfg,figCfg,evaluators,handles) ...
    buildSweeps(cfg, figCfg, evaluators, handles, modelFns, mValues, lambdaVals);
end

function lambdaVals = applyLambdaExclusions(lambdaVals, ranges)
if nargin < 2 || isempty(ranges)
    return
end
if ~isnumeric(ranges) || size(ranges,2) ~= 2
    return
end
mask = true(size(lambdaVals));
for k = 1:size(ranges,1)
    lo = min(ranges(k,1), ranges(k,2));
    hi = max(ranges(k,1), ranges(k,2));
    mask = mask & ~(lambdaVals >= lo & lambdaVals <= hi);
end
lambdaVals = lambdaVals(mask);
end

function handles = buildHandles(cfg, fns)
handles = struct();
handles.ode = @(x,y) fns.ode(x,y, cfg.p, cfg.n);
handles.bc  = @(ya,yb) fns.bc(ya,yb, cfg.p, cfg.n);
end

function metricEvaluators = buildMetricEvaluators(cfg, fns)
metricEvaluators = struct( ...
    'Cf_star', @(sol) fns.localSkinFriction(sol, cfg.p, cfg.n), ...
    'Nu_star', @(sol) fns.nusseltNumber(sol, cfg.p, cfg.n));
end

function sweeps = buildSweeps(cfg, figureConfigs, metricEvaluators, handles, fns, mValues, lambdaVals)
if nargin < 4 || isempty(handles)
    handles = buildHandles(cfg, fns);
end %#ok<NASGU>

sweeps = struct('name', {}, 'options', {});

mSweep = struct();
mSweep.name = 'M';
mSweep.options = struct();
mSweep.options.sweep = struct( ...
    'name', 'M', ...
    'values', mValues, ...
    'odeFactory', @(val) parameterizedOde(cfg, 'M', val, fns.ode), ...
    'bcFactory', @(val) parameterizedBc(cfg, 'M', val, fns.bc), ...
    'labelFcn', @(val) sprintf('M=%.2f', val));
mSweep.options.figureConfigs = figureConfigs(1:2);
mSweep.options.metricEvaluators = metricEvaluators;
mSweep.options.domainLabel = 'M';
mSweep.optionsBuilder = [];
mSweep.groupParamName = '';
mSweep.groupValues = [];
mSweep.groupLabelFcn = [];

lambdaSweep = struct();
lambdaSweep.name = '\lambda';
lambdaSweep.optionsBuilder = @(cfgLocal, figCfg, metrics) buildLambdaOptions(cfgLocal, figCfg, metrics, fns, lambdaVals);
lambdaSweep.options = lambdaSweep.optionsBuilder(cfg, figureConfigs(3:4), metricEvaluators);
lambdaSweep.options.domainLabel = '\lambda';
lambdaSweep.groupParamName = 'M';
lambdaSweep.groupValues = mValues;
lambdaSweep.groupLabelFcn = @(val) sprintf('M=%.2f', val);

sweeps = [mSweep, lambdaSweep];
end

function opts = buildLambdaOptions(cfg, figureConfigs, metricEvaluators, fns, lambdaVals)
sweepStruct = struct( ...
    'name', '\lambda', ...
    'values', lambdaVals, ...
    'startValue', lambdaVals(1), ...
    'odeFactory', @(val) parameterizedOde(cfg, 'lambda', val, fns.ode), ...
    'bcFactory', @(val) parameterizedBc(cfg, 'lambda', val, fns.bc), ...
    'labelFcn', @(val) sprintf('\\lambda=%.2f', val), ...
    'probeSides', true, ...
    'maxConsecutiveFails', 1);
opts = struct();
opts.sweep = sweepStruct;
opts.refine = struct('onFail', true, 'numSubdiv', 5);
opts.figureConfigs = figureConfigs;
opts.metricEvaluators = metricEvaluators;
end

function configs = buildFigureConfigs(cfg, metricEvaluators)
cfEvaluator = pickMetric(metricEvaluators, 'Cf_star');
nuEvaluator = pickMetric(metricEvaluators, 'Nu_star');
plotOutliers = false;
if isfield(cfg,'plotOutliers') && ~isempty(cfg.plotOutliers)
    plotOutliers = logical(cfg.plotOutliers);
end
asymptoteCfg = struct( ...
    'enabled', true, ...
    'slopeFactor', 10, ...
    'yPercentile', 95, ...
    'yMadFactor', 8, ...
    'yJumpFactor', 10, ...
    'fadeFactor', 0.5, ...
    'breakLine', true);

configs = [
    struct('mode','profile','figureId',1,'componentIdx',2,'domainWindow',[0,10],'codomainWindow',[],'numSamples',800, ...
    'xlabel','\eta','ylabel','f`(\eta)','title','Figure 1: Velocity Profile', 'showOutliers', plotOutliers, ...
    'metricField','Cf_star','metricFn',[],'branchIdx',1:2,'useSharedColor',false, 'asymptote', struct('enabled', false));
    struct('mode','profile','figureId',2,'componentIdx',4,'domainWindow',[0,10],'codomainWindow',[],'numSamples',801, ...
    'xlabel','\eta','ylabel','\theta(\eta)','title','Figure 2: Temperature Profile', 'showOutliers', plotOutliers, ...
    'metricField','Nu_star','metricFn',[],'branchIdx',1:2,'useSharedColor',false, 'asymptote', struct('enabled', false));
    struct('mode','metric','figureId',3,'componentIdx',[],'domainWindow',[],'codomainWindow',[], 'numSamples',800, ...
    'xlabel','\lambda','ylabel','C_f R_e^{1/2}','title','Figure 3: Skin friction vs \lambda', 'showOutliers', plotOutliers, ...
    'metricField','Cf_star','metricFn',cfEvaluator,'branchIdx',1:2,'useSharedColor',true, 'asymptote', struct('enabled', false));
    struct('mode','metric','figureId',4,'componentIdx',[],'domainWindow',[],'codomainWindow',[],'numSamples',800, ...
    'xlabel','\lambda','ylabel','Nu R_e^{-1/2}','title','Figure 4: Nusselt vs \lambda', 'showOutliers', plotOutliers, ...
    'metricField','Nu_star','metricFn',nuEvaluator,'branchIdx',1:2,'useSharedColor',true, 'asymptote', asymptoteCfg)
    ];
end

function fn = pickMetric(metricEvaluators, fieldName)
fn = [];
if isempty(metricEvaluators)
    return
end
if isfield(metricEvaluators, fieldName) && ~isempty(metricEvaluators.(fieldName))
    fn = metricEvaluators.(fieldName);
end
end

function handle = parameterizedOde(cfg, fieldName, value, odeFn)
pTmp = cfg.p;
pTmp.(fieldName) = value;
handle = @(x,y) odeFn(x,y,pTmp,cfg.n);
end

function handle = parameterizedBc(cfg, fieldName, value, bcFn)
pTmp = cfg.p;
pTmp.(fieldName) = value;
handle = @(ya,yb) bcFn(ya,yb,pTmp,cfg.n);
end
