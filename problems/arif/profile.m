function profileSpec = profile()
%PROFILE Sweep definitions and evaluators for Natasha.

modelFns = model();

mValues = [0.1, 0.2, 0.3];
lambdaVals = linspace(-1.8, 1.5, 30);

profileSpec = struct();
profileSpec.createBaseHandles = @(cfg) buildHandles(cfg, modelFns);
profileSpec.metricEvaluators  = @(cfg) buildMetricEvaluators(cfg, modelFns);
profileSpec.figureConfigs     = @(cfg, metrics) buildFigureConfigs(cfg, metrics);
profileSpec.sweeps            = @(cfg,figCfg,evaluators,handles) ...
    buildSweeps(cfg, figCfg, evaluators, handles, modelFns, mValues, lambdaVals);
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
    'labelFcn', @(val) sprintf('M=%.1f', val));
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
lambdaSweep.groupLabelFcn = @(val) sprintf('M=%.1f', val);

sweeps = [mSweep, lambdaSweep];
end

function opts = buildLambdaOptions(cfg, figureConfigs, metricEvaluators, fns, lambdaVals)
sweepStruct = struct( ...
    'name', '\lambda', ...
    'values', lambdaVals, ...
    'startValue', lambdaVals(1), ...
    'odeFactory', @(val) parameterizedOde(cfg, 'lambda', val, fns.ode), ...
    'bcFactory', @(val) parameterizedBc(cfg, 'lambda', val, fns.bc), ...
    'labelFcn', @(val) sprintf('\\lambda=%.2f', val));
opts = struct();
opts.sweep = sweepStruct;
opts.refine = struct('onFail', true, 'numSubdiv', 5);
opts.figureConfigs = figureConfigs;
opts.metricEvaluators = metricEvaluators;
end

function configs = buildFigureConfigs(cfg, metricEvaluators)
cfEvaluator = pickMetric(metricEvaluators, 'Cf_star');
nuEvaluator = pickMetric(metricEvaluators, 'Nu_star');

configs = [
    struct('mode','profile','figureId',1,'componentIdx',2,'domainWindow',[0,3.5],'codomainWindow',[],'numSamples',800, ...
    'xlabel','\eta','ylabel','f`(\eta)','title','Figure 1: Velocity Profile', ...
    'metricField','Cf_star','metricFn',[],'branchIdx',1:2,'useSharedColor',false);
    struct('mode','profile','figureId',2,'componentIdx',4,'domainWindow',[0,2],'codomainWindow',[],'numSamples',801, ...
    'xlabel','\eta','ylabel','\t(\eta)','title','Figure 2: Temperature Profile', ...
    'metricField','Nu_star','metricFn',[],'branchIdx',1:2,'useSharedColor',false);
    struct('mode','metric','figureId',3,'componentIdx',[],'domainWindow',[],'codomainWindow',[], 'numSamples',800, ...
    'xlabel','\lambda','ylabel','C_f R_e^{1/2}','title','Figure 3: Skin friction vs \lambda', ...
    'metricField','Cf_star','metricFn',cfEvaluator,'branchIdx',1:2,'useSharedColor',true);
    struct('mode','metric','figureId',4,'componentIdx',[],'domainWindow',[],'codomainWindow',[],'numSamples',800, ...
    'xlabel','\lambda','ylabel','Nu R_e^{1/2}','title','Figure 4: Nusselt vs \lambda', ...
    'metricField','Nu_star','metricFn',nuEvaluator,'branchIdx',1:2,'useSharedColor',true)
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
