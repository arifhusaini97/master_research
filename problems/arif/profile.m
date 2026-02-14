function profileSpec = profile()
%PROFILE Sweep definitions and evaluators for Arif.

modelFns = model();
vals = values();

%No 2nd solution for M=0
mValues = vals.sweep.mValues;
excludeRanges = vals.sweep.excludeLambdaRanges;
lambdaVals = [];
if isfield(vals.sweep,'lambdaVals') && ~isempty(vals.sweep.lambdaVals)
    lambdaVals = applyLambdaExclusions(vals.sweep.lambdaVals, excludeRanges);
end
smartCfg = [];
if isfield(vals.sweep,'smartLinspaceLambdaVals') && ~isempty(vals.sweep.smartLinspaceLambdaVals) ...
        && isstruct(vals.sweep.smartLinspaceLambdaVals)
    smartCfg = vals.sweep.smartLinspaceLambdaVals;
end
% Default if no values.sweep.primary and values.sweep.secondary in values.m
primarySpec = resolveSweepSpec(vals.sweep, 'primary', struct( ...
    'paramScope', 'p', ...
    'paramName', 'M', ...
    'name', 'M', ...
    'values', mValues, ...
    'labelFcn', @(val) sprintf('M=%.2f', val)));
secondarySpec = resolveSweepSpec(vals.sweep, 'secondary', struct( ...
    'paramScope', 'p', ...
    'paramName', 'lambda', ...
    'name', '\lambda', ...
    'values', lambdaVals, ...
    'labelFcn', @(val) sprintf('\\lambda=%.2f', val), ...
    'probeSides', true, ...
    'maxConsecutiveFails', 1, ...
    'refine', struct('onFail', true, 'numSubdiv', 5), ...
    'useSmartLinspace', true, ...
    'excludeRanges', excludeRanges, ...
    'group', struct( ...
    'paramScope', 'p', ...
    'paramName', 'M', ...
    'values', mValues, ...
    'labelFcn', @(val) sprintf('M=%.2f', val))));
%%%%
profileSpec = struct();
profileSpec.createBaseHandles = @(cfg) buildHandles(cfg, modelFns);
profileSpec.metricEvaluators  = @(cfg) buildMetricEvaluators(cfg, modelFns);
profileSpec.figureConfigs     = @(cfg, metrics) buildFigureConfigs(cfg, metrics);
profileSpec.sweeps            = @(cfg,figCfg,evaluators,handles) ...
    buildSweeps(cfg, figCfg, evaluators, handles, modelFns, primarySpec, secondarySpec, smartCfg, excludeRanges);
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
    'Cf_star', @(sol, varargin) evaluateMetricWithSweep(sol, cfg, fns, varargin, 'Cf_star'), ...
    'Nu_star', @(sol, varargin) evaluateMetricWithSweep(sol, cfg, fns, varargin, 'Nu_star'));
end

function sweeps = buildSweeps(cfg, figureConfigs, metricEvaluators, handles, fns, primarySpec, secondarySpec, smartCfg, excludeRanges)
if nargin < 4 || isempty(handles)
    handles = buildHandles(cfg, fns);
end %#ok<NASGU>

sweeps = struct('name', {}, 'options', {});

primarySweep = buildSweepFromSpec(cfg, figureConfigs(1:2), metricEvaluators, fns, primarySpec, smartCfg, excludeRanges);
secondarySweep = buildSweepFromSpec(cfg, figureConfigs(3:4), metricEvaluators, fns, secondarySpec, smartCfg, excludeRanges);

sweeps = [primarySweep, secondarySweep];
end

function sigma = getOutlierSigma(cfg)
sigma = 3;
if isfield(cfg,'outlierSigma') && ~isempty(cfg.outlierSigma)
    candidate = cfg.outlierSigma;
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && candidate > 0
        sigma = candidate;
    end
end
end

function opts = buildSweepOptions(cfg, figureConfigs, metricEvaluators, fns, spec, smartCfg, excludeRanges)
valuesToUse = resolveSweepValues(cfg, spec, smartCfg, excludeRanges);
if isempty(valuesToUse)
    error('profile:MissingSweepVals', 'No sweep values configured for %s.', resolveParamLabel(spec));
end
labelFcn = resolveLabelFcn(spec);
sweepStruct = struct( ...
    'name', resolveSweepName(spec), ...
    'values', valuesToUse, ...
    'startValue', valuesToUse(1), ...
    'paramScope', spec.paramScope, ...
    'paramName', spec.paramName, ...
    'odeFactory', @(val) parameterizedOde(cfg, spec, val, fns.ode), ...
    'bcFactory', @(val) parameterizedBc(cfg, spec, val, fns.bc), ...
    'labelFcn', labelFcn);
if isfield(spec,'probeSides')
    sweepStruct.probeSides = spec.probeSides;
end
if isfield(spec,'maxConsecutiveFails')
    sweepStruct.maxConsecutiveFails = spec.maxConsecutiveFails;
end
sweepStruct.smartRateEnabled = isfield(cfg,'smartLinspace') && cfg.smartLinspace;
sweepStruct.smartRateConfig = smartCfg;
opts = struct();
opts.sweep = sweepStruct;
if isfield(spec,'refine') && isstruct(spec.refine)
    opts.refine = spec.refine;
end
opts.figureConfigs = figureConfigs;
opts.metricEvaluators = metricEvaluators;
opts.outlierSigma = getOutlierSigma(cfg);
opts.domainLabel = resolveSweepLabel(spec);
end

function configs = buildFigureConfigs(cfg, metricEvaluators)
cfEvaluator = pickMetric(metricEvaluators, 'Cf_star');
nuEvaluator = pickMetric(metricEvaluators, 'Nu_star');
plotOutliers = false;
if isfield(cfg,'plotOutliers') && ~isempty(cfg.plotOutliers)
    plotOutliers = logical(cfg.plotOutliers);
end
outlierFilter = 'global';
if isfield(cfg,'outlierFilterMode') && ~isempty(cfg.outlierFilterMode)
    outlierFilter = cfg.outlierFilterMode;
end
plotZeroDeviation = true;
if isfield(cfg,'plotZeroDeviationPoints') && ~isempty(cfg.plotZeroDeviationPoints)
    plotZeroDeviation = logical(cfg.plotZeroDeviationPoints);
end
zeroDeviationTol = 1e-6;
if isfield(cfg,'zeroDeviationTolerance') && ~isempty(cfg.zeroDeviationTolerance)
    zeroDeviationTol = cfg.zeroDeviationTolerance;
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
    'metricField','Cf_star','metricFn',[],'branchIdx',1:2,'useSharedColor',false, 'outlierFilter', outlierFilter, ...
    'plotZeroDeviation', plotZeroDeviation, 'zeroDeviationTolerance', zeroDeviationTol, ...
    'asymptote', struct('enabled', false));
    struct('mode','profile','figureId',2,'componentIdx',4,'domainWindow',[0,10],'codomainWindow',[],'numSamples',801, ...
    'xlabel','\eta','ylabel','\theta(\eta)','title','Figure 2: Temperature Profile', 'showOutliers', plotOutliers, ...
    'metricField','Nu_star','metricFn',[],'branchIdx',1:2,'useSharedColor',false, 'outlierFilter', outlierFilter, ...
    'plotZeroDeviation', plotZeroDeviation, 'zeroDeviationTolerance', zeroDeviationTol, ...
    'asymptote', struct('enabled', false));
    struct('mode','metric','figureId',3,'componentIdx',[],'domainWindow',[],'codomainWindow',[], 'numSamples',800, ...
    'xlabel','\lambda','ylabel','C_f R_e^{1/2}','title','Figure 3: Skin friction vs Bouyancy parameter', 'showOutliers', plotOutliers, ...
    'metricField','Cf_star','metricFn',cfEvaluator,'branchIdx',1:2,'useSharedColor',true, 'outlierFilter', outlierFilter, ...
    'plotZeroDeviation', plotZeroDeviation, 'zeroDeviationTolerance', zeroDeviationTol, ...
    'asymptote', struct('enabled', false));
    struct('mode','metric','figureId',4,'componentIdx',[],'domainWindow',[],'codomainWindow',[],'numSamples',800, ...
    'xlabel','\lambda','ylabel','Nu R_e^{-1/2}','title','Figure 4: Nusselt vs Bouyancy parameter', 'showOutliers', plotOutliers, ...
    'metricField','Nu_star','metricFn',nuEvaluator,'branchIdx',1:2,'useSharedColor',true, 'outlierFilter', outlierFilter, ...
    'plotZeroDeviation', plotZeroDeviation, 'zeroDeviationTolerance', zeroDeviationTol, ...
    'asymptote', asymptoteCfg)
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

function handle = parameterizedOde(cfg, spec, value, odeFn)
[pTmp, ~, nTmp] = applyParamOverride(cfg, spec, value);
handle = @(x,y) odeFn(x,y,pTmp,nTmp);
end

function handle = parameterizedBc(cfg, spec, value, bcFn)
[pTmp, ~, nTmp] = applyParamOverride(cfg, spec, value);
handle = @(ya,yb) bcFn(ya,yb,pTmp,nTmp);
end

function specOut = resolveSweepSpec(sweepCfg, fieldName, defaults)
specOut = defaults;
if nargin < 1 || isempty(sweepCfg)
    return
end
if ~isfield(sweepCfg, fieldName)
    return
end
candidate = sweepCfg.(fieldName);
if ~isstruct(candidate)
    return
end
names = fieldnames(candidate);
for k = 1:numel(names)
    specOut.(names{k}) = candidate.(names{k});
end
end

function sweepDef = buildSweepFromSpec(cfg, figCfg, metrics, fns, spec, smartCfg, excludeRanges)
if nargin < 2 || isempty(figCfg)
    figCfg = [];
end
sweepDef = struct();
sweepDef.name = resolveSweepName(spec);
if isfield(spec,'group') && isstruct(spec.group) && isfield(spec.group,'values') && ~isempty(spec.group.values)
    sweepDef.optionsBuilder = @(cfgLocal, figCfgLocal, metricsLocal) buildSweepOptions(cfgLocal, figCfgLocal, metricsLocal, fns, spec, smartCfg, excludeRanges);
    sweepDef.options = buildSweepOptions(cfg, figCfg, metrics, fns, spec, smartCfg, excludeRanges);
    sweepDef.groupParamName = spec.group.paramName;
    sweepDef.groupParamScope = resolveParamScope(spec.group);
    sweepDef.groupValues = spec.group.values;
    sweepDef.groupLabelFcn = resolveLabelFcn(spec.group);
else
    sweepDef.options = buildSweepOptions(cfg, figCfg, metrics, fns, spec, smartCfg, excludeRanges);
    sweepDef.optionsBuilder = [];
    sweepDef.groupParamName = '';
    sweepDef.groupParamScope = '';
    sweepDef.groupValues = [];
    sweepDef.groupLabelFcn = [];
end
end

function valuesToUse = resolveSweepValues(cfg, spec, smartCfg, excludeRanges)
valuesToUse = [];
if isfield(spec,'values') && ~isempty(spec.values)
    valuesToUse = spec.values;
end
useSmart = isfield(spec,'useSmartLinspace') && spec.useSmartLinspace;
if useSmart && isfield(cfg,'smartLinspace') && cfg.smartLinspace ...
        && isstruct(smartCfg) && isfield(smartCfg,'min') && isfield(smartCfg,'max') && isfield(smartCfg,'count')
    valuesToUse = linspace(smartCfg.min, smartCfg.max, smartCfg.count);
end
exclusions = excludeRanges;
if isfield(spec,'excludeRanges') && ~isempty(spec.excludeRanges)
    exclusions = spec.excludeRanges;
end
if ~isempty(valuesToUse) && ~isempty(exclusions)
    valuesToUse = applyLambdaExclusions(valuesToUse, exclusions);
end
end

function name = resolveSweepName(spec)
name = resolveParamLabel(spec);
if isfield(spec,'name') && strlength(string(spec.name)) > 0
    name = spec.name;
end
end

function label = resolveSweepLabel(spec)
label = resolveSweepName(spec);
if isfield(spec,'domainLabel') && strlength(string(spec.domainLabel)) > 0
    label = spec.domainLabel;
end
end

function label = resolveParamLabel(spec)
label = '';
if isfield(spec,'paramName') && strlength(string(spec.paramName)) > 0
    label = spec.paramName;
end
end

function scope = resolveParamScope(spec)
scope = 'p';
if isfield(spec,'paramScope') && strlength(string(spec.paramScope)) > 0
    scope = char(lower(string(spec.paramScope)));
end
end

function labelFcn = resolveLabelFcn(spec)
if isfield(spec,'labelFcn') && ~isempty(spec.labelFcn)
    labelFcn = spec.labelFcn;
    return
end
paramLabel = resolveSweepLabel(spec);
labelFcn = @(val) sprintf('%s=%s', char(paramLabel), formatNumericToken(val));
end

function [pTmp, mTmp, nTmp] = applyParamOverride(cfg, spec, value)
pTmp = cfg.p;
mTmp = cfg.m;
nTmp = cfg.n;
scope = resolveParamScope(spec);
name = spec.paramName;
if strcmp(scope, 'm')
    mTmp.(name) = value;
    if isfield(cfg,'deriveNFromM') && isa(cfg.deriveNFromM, 'function_handle')
        nTmp = cfg.deriveNFromM(mTmp);
    end
else
    pTmp.(name) = value;
end
end

function value = evaluateMetricWithSweep(sol, cfg, fns, contextArgs, metricName)
ctx = [];
if nargin >= 4 && ~isempty(contextArgs)
    ctx = contextArgs{1};
end
[pTmp, ~, nTmp] = resolveMetricParams(cfg, ctx);
if strcmpi(metricName, 'Nu_star')
    value = fns.nusseltNumber(sol, pTmp, nTmp);
else
    value = fns.localSkinFriction(sol, pTmp, nTmp);
end
end

function [pTmp, mTmp, nTmp] = resolveMetricParams(cfg, ctx)
pTmp = cfg.p;
mTmp = cfg.m;
nTmp = cfg.n;
if isempty(ctx) || ~isstruct(ctx)
    return
end
if ~isfield(ctx,'paramName') || ~isfield(ctx,'paramScope') || ~isfield(ctx,'paramValue')
    return
end
spec = struct('paramName', ctx.paramName, 'paramScope', ctx.paramScope);
[pTmp, mTmp, nTmp] = applyParamOverride(cfg, spec, ctx.paramValue);
end

function token = formatNumericToken(value)
token = numericFormat('token', value);
end
