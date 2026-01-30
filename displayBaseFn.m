% displayBaseFn.m
function [results, sweepMeta, figurePayloads] = displayBaseFn(odeFun, bcFun, guesses, domainMinList, domainMaxList, gridSizeList, methodCfg, varargin)
% displayBaseFn  Solve and plot multiple BVP branches with optional sweeps.
% odeFun, bcFun : baseline function handles with params captured
% guesses       : cell array of @(x) initial-guess functions, one per branch
% domainMinList : scalar start of domain (shared) or vector per-branch
% domainMaxList : scalar or vector of per-branch domain maxima
% gridSizeList  : scalar or vector of per-branch grid sizes (steps)
% varargin      : optional logical doPlot flag or struct of options

if nargin < 7 || isempty(methodCfg)
    methodCfg = struct();
end
methodCfg = normalizeMethodSpec(methodCfg, domainMinList, domainMaxList);

opts = parseDisplayOptions(varargin);
doPlot   = opts.doPlot;
sweepCfg = opts.sweep;
figConfigs = opts.figureConfigs;
metricEvaluators = opts.metricEvaluators;

[values, odeHandles, bcHandles, sweepLabels, sweepName, sweepCfgOut] = ...
    buildSweepHandles(sweepCfg, odeFun, bcFun);
if isfield(sweepCfgOut,'paramName') && strlength(string(sweepCfgOut.paramName)) > 0
    methodCfg.paramName = sweepCfgOut.paramName;
end
if isfield(sweepCfgOut,'paramScope') && strlength(string(sweepCfgOut.paramScope)) > 0
    methodCfg.paramScope = sweepCfgOut.paramScope;
end

nBranches  = numel(guesses);
blankTemplate = blankResultStruct();
numSweeps  = numel(odeHandles);
results    = repmat(blankTemplate, nBranches, numSweeps);
logPlannedPointCount(sweepName, values, nBranches);

prevSolutions = cell(1, nBranches);
solutionCache = cell(nBranches, numSweeps);
processed = false(1, numSweeps);
hasSuccessfulPoint = false;
probeCfg = resolveProbeConfig(sweepCfgOut);
probeState = initProbeState(probeCfg);
rightToLeftFailCount = 0;
attemptLog = emptyAttemptLog();
lineBlockIdx = beginLineLog();

while true
    if isempty(values) || isempty(processed)
        break;
    end
    [sweepIdx, probeState] = selectNextSweepIndex(values, processed, sweepCfgOut, probeState);
    if isempty(sweepIdx)
        break;
    end

    currentValue = pickValueSafe(values, sweepIdx);
    if isempty(sweepLabels)
        currentLabel = '';
    else
        currentLabel = sweepLabels{min(sweepIdx, numel(sweepLabels))};
    end

    [results, solutionCache, prevSolutions, branchSuccess, branchDiagnostics] = ...
        solveSweepPoint(sweepIdx, values, odeHandles, bcHandles, guesses, domainMinList, domainMaxList, gridSizeList, ...
        methodCfg, results, solutionCache, prevSolutions, odeFun, bcFun, metricEvaluators, opts.domainLabel);

    sweepStep = computeSweepStep(values, sweepIdx);
    branchDiagnostics = assignSweepStep(branchDiagnostics, sweepStep);
    attemptLog(end+1) = buildAttemptLogEntry(currentValue, currentLabel, branchDiagnostics); %#ok<AGROW>

    successThisPoint = any(branchSuccess);
    processed(sweepIdx) = true;
    if successThisPoint
        hasSuccessfulPoint = true;
        if strcmp(probeState.lastProbeSide, 'right')
            probeState.rightSuccess = true;
        end
    else
    end
    probeState = updateProbeStateAfterAttempt(probeState, values, processed);
    if strcmp(probeState.mode, 'right_to_left')
        if successThisPoint
            rightToLeftFailCount = 0;
        else
            rightToLeftFailCount = rightToLeftFailCount + 1;
            if rightToLeftFailCount >= probeCfg.maxConsecutiveFails && probeCfg.maxConsecutiveFails > 0
                emitRunLog('Right-to-left sweep: %d consecutive fails; skipping remaining sweep values.\n', rightToLeftFailCount);
                attemptLog = appendSkipAttempts(attemptLog, values, sweepLabels, processed, nBranches, ...
                    'Right-to-left sweep skipped remaining values after consecutive failures.');
                break;
            end
        end
    end

    if ~successThisPoint
        if shouldRefine(sweepCfgOut) && hasSuccessfulPoint && ~strcmp(probeState.mode, 'right_to_left')
            [values, odeHandles, bcHandles, sweepLabels, results, solutionCache, processed, prevSolutions, newAttempts] = ...
                handleFailedSweepPoint(sweepIdx, values, odeHandles, bcHandles, sweepLabels, results, solutionCache, ...
                processed, prevSolutions, guesses, domainMinList, domainMaxList, gridSizeList, ...
                methodCfg, sweepCfgOut, odeFun, bcFun, metricEvaluators, opts.domainLabel);
            if ~isempty(newAttempts)
                attemptLog = [attemptLog, newAttempts]; %#ok<AGROW>
            end
            if isempty(values)
                break;
            end
        else
            [values, odeHandles, bcHandles, sweepLabels, results, solutionCache, processed] = ...
                removeValueAtIndex(values, odeHandles, bcHandles, sweepLabels, results, solutionCache, processed, sweepIdx);
            if isempty(values)
                break;
            end
        end
    end
end

attemptLog = computeAttemptRates(attemptLog);
 [values, sweepLabels, results, attemptLog] = applySmartRateFilter(values, sweepLabels, results, attemptLog, sweepCfgOut);
failRanges = computeContinuousFailRanges(attemptLog, 'fail');
successRanges = computeContinuousFailRanges(attemptLog, 'success');
successOutliers = detectSuccessOutliers(attemptLog, opts.outlierSigma);
metricFallbacks = buildMetricFallbacks(attemptLog);
logContinuousRangeSummary(failRanges, sweepName, opts.domainLabel, 'fail');
logContinuousRangeSummary(successRanges, sweepName, opts.domainLabel, 'success');
logOutlierSummary(successOutliers, sweepName, opts.domainLabel);
endLineLog(lineBlockIdx);

numSweeps = size(results, 2);
colors     = lines(max(numSweeps, 1));
lineStyles = {'-','--',':','-.'};

if numSweeps == 1
    results = results(:,1);
end

sweepMeta = struct( ...
    'sweepValues', values, ...
    'sweepLabels', {sweepLabels}, ...
    'sweepName', sweepName, ...
    'figureConfigs', figConfigs, ...
    'rawSweepConfig', sweepCfgOut, ...
    'attemptLog', attemptLog, ...
    'continuousFailRanges', failRanges, ...
    'continuousSuccessRanges', successRanges, ...
    'successOutliers', successOutliers, ...
    'domainLabel', opts.domainLabel);

figurePayloads = {};
if doPlot || nargout > 2
figurePayloads = renderFigures(results, values, sweepLabels, colors, lineStyles, figConfigs, attemptLog, failRanges, successRanges, successOutliers, metricFallbacks, doPlot);
end
end

% --- helper: pick ith value from scalar or vector (fallback to last) ---
function v = pick(list, i)
if isscalar(list)
    v = list;
else
    v = list(min(i, numel(list)));
end
end

function opts = parseDisplayOptions(args)
opts = struct('doPlot', true, 'sweep', [], 'figureConfigs', [], 'metricEvaluators', struct(), ...
    'domainLabel', 'x', 'outlierSigma', 3);
if isempty(args)
    return
end
for k = 1:numel(args)
    candidate = args{k};
    if islogical(candidate)
        opts.doPlot = candidate;
    elseif isstruct(candidate)
        if isfield(candidate,'doPlot')
            opts.doPlot = candidate.doPlot;
        end
        if isfield(candidate,'sweep')
            opts.sweep = candidate.sweep;
        end
        if isfield(candidate,'figureConfigs')
            opts.figureConfigs = candidate.figureConfigs;
        end
        if isfield(candidate,'metricEvaluators')
            opts.metricEvaluators = candidate.metricEvaluators;
        end
        if isfield(candidate,'refine') && ~isempty(candidate.refine)
            opts.sweep.refine = candidate.refine;
        end
        if isfield(candidate,'domainLabel') && ~isempty(candidate.domainLabel)
            opts.domainLabel = candidate.domainLabel;
        end
        if isfield(candidate,'outlierSigma') && ~isempty(candidate.outlierSigma)
            opts.outlierSigma = candidate.outlierSigma;
        end
    end
end
if isempty(opts.figureConfigs)
    opts.figureConfigs = struct('figureId',1,'componentIdx',2,'domainWindow',[0,20],'codomainWindow',[], ...
        'numSamples',800,'xlabel','x','ylabel','codomain','title','Figure 1');
end
end

function [values, odeHandles, bcHandles, labels, name, cfgOut] = buildSweepHandles(cfg, baseOde, baseBc)
if isempty(cfg)
    values     = NaN;
    odeHandles = {baseOde};
    bcHandles  = {baseBc};
    labels     = {''};
    name       = '';
    cfgOut     = struct();
    return
end

if isfield(cfg,'values')
    values = cfg.values;
elseif isfield(cfg,'range') && isfield(cfg,'step')
    values = cfg.range(1):cfg.step:cfg.range(2);
else
    values = NaN;
end
cfgOut = cfg;

name = '';
if isfield(cfg,'name')
    name = cfg.name;
end

if isfield(cfg,'labelFcn')
    labelFcn = cfg.labelFcn;
else
    labelFcn = @(val) sprintf('%s=%s', name, formatNumericToken(val));
end

odeHandles = cell(1, numel(values));
bcHandles  = cell(1, numel(values));
labels     = cell(1, numel(values));

for k = 1:numel(values)
    val = values(k);
    if isfield(cfg,'odeFactory')
        odeHandles{k} = cfg.odeFactory(val);
    else
        odeHandles{k} = baseOde;
    end
    if isfield(cfg,'bcFactory')
        bcHandles{k} = cfg.bcFactory(val);
    else
        bcHandles{k} = baseBc;
    end
    labels{k} = labelFcn(val);
end
end

function [sol, diagnostics] = solveBranch(odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel)
solverId = lower(strtrim(string(methodCfg.solver)));
solverFn = resolveSolverFunction(solverId);
utils = buildSolverUtilities();
[sol, diagnostics] = solverFn(odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel, utils);
end

function solverFn = resolveSolverFunction(solverId)
% Objective: Provide a single lookup for solver implementations.
% Purpose : Decouple display logic from per-solver plumbing.
% SWOT - S: central registry. W: requires valid solverId. O: easy new solvers. T: missing files break resolves.
%
persistent registry methodsPath
if isempty(registry)
    if isempty(methodsPath)
        methodsPath = fullfile(fileparts(mfilename('fullpath')), 'methods');
        if exist(methodsPath, 'dir')
            addpath(methodsPath);
        end
    end
    registry = containers.Map( ...
        {'bvp4c','collocation4','bvp5c','collocation5','bdf2'}, ...
        {@bvp4cSolver,@bvp4cSolver,@bvp5cSolver,@bvp5cSolver,@bdf2Solver});
end
if ~registry.isKey(solverId)
    error('displayBaseFn:unknownMethod','Unknown solver "%s".', solverId);
end
solverFn = registry(solverId);
end

function utils = buildSolverUtilities()
utils = struct();
utils.computeInitialGuessError = @computeGuessErrorWrapper;
utils.averageMeshSpacing = @averageSpacingWrapper;
utils.extractMaxResidual = @extractMaxResidualWrapper;
utils.extractIterationCount = @extractIterationWrapper;
utils.formatSolverConsoleLog = @formatConsoleWrapper;
utils.buildMethodConsoleLog = @buildConsoleBodyWrapper;
utils.enforceTimeLimit = @enforceTimeLimitWrapper;
utils.midpoints = @midpointsWrapper;
utils.initBranchDiagnostics = @initDiagnosticsWrapper;
utils.getfieldWithDefault = @getfieldWithDefaultWrapper;
utils.beginLineLog = @beginLineLogWrapper;
utils.endLineLog = @endLineLogWrapper;
utils.emitRunLog = @emitRunLog;
utils.sampleInitialGuess = @sampleInitialGuessWrapper;

    function err = computeGuessErrorWrapper(sol, guess)
        err = computeInitialGuessError(sol, guess);
    end
    function avg = averageSpacingWrapper(mesh)
        avg = averageMeshSpacing(mesh);
    end
    function res = extractMaxResidualWrapper(sol)
        res = extractMaxResidual(sol);
    end
    function iters = extractIterationWrapper(sol)
        iters = extractIterationCount(sol);
    end
    function txt = formatConsoleWrapper(rawLog, attemptNum, domainMin, domainMax, domainLabel, varargin)
        meta = [];
        if ~isempty(varargin)
            meta = varargin{1};
        end
        txt = formatSolverConsoleLog(rawLog, attemptNum, domainMin, domainMax, domainLabel, meta);
    end
    function body = buildConsoleBodyWrapper(meshPts, residual, odeCalls, bcCalls)
        body = buildMethodConsoleLog(meshPts, residual, odeCalls, bcCalls);
    end
    function out = enforceTimeLimitWrapper(fun, startClock, timeoutSeconds, varargin)
        out = enforceTimeLimit(fun, startClock, timeoutSeconds, varargin{:});
    end
    function mids = midpointsWrapper(mesh)
        mids = midpoints(mesh);
    end
    function diag = initDiagnosticsWrapper(stepSize)
        diag = initBranchDiagnostics(stepSize);
    end
    function val = getfieldWithDefaultWrapper(s, fieldName, defaultVal)
        val = getfieldWithDefault(s, fieldName, defaultVal);
    end
    function idx = beginLineLogWrapper()
        idx = beginLineLog();
    end
    function endLineLogWrapper(idx)
        endLineLog(idx);
    end
    function vals = sampleInitialGuessWrapper(guessFcn, mesh)
        vals = evaluateGuessOnMesh(guessFcn, mesh);
    end
end
function sol = selectNearestSolution(cacheRow, values, targetIdx)
sol = [];
if isempty(values) || isnan(values(1)) || targetIdx > numel(values)
    return
end
diffs = inf(size(values));
for k = 1:numel(cacheRow)
    if ~isempty(cacheRow{k})
        diffs(k) = abs(values(k) - values(targetIdx));
    end
end
[~, idx] = min(diffs);
if isfinite(diffs(idx))
    sol = cacheRow{idx};
end
end

function order = determineSolveOrder(values, cfg)
n = numel(values);
order = 1:n;
if n == 0
    order = [];
    return
end
if nargin < 2 || isempty(cfg) || ~isfield(cfg,'startValue') || isempty(cfg.startValue) || isnan(values(1))
    return
end
[~,startIdx] = min(abs(values - cfg.startValue));
order = startIdx;
left = startIdx - 1;
right = startIdx + 1;
while left >= 1 || right <= n
    if right <= n
        order(end+1) = right; %#ok<AGROW>
        right = right + 1;
    end
    if left >= 1
        order(end+1) = left; %#ok<AGROW>
        left = left - 1;
    end
end
end

function order = remainingSolveOrder(values, processed, cfg)
fullOrder = determineSolveOrder(values, cfg);
mask = ~processed(fullOrder);
order = fullOrder(mask);
end

function stepVal = computeSweepStep(values, idx)
stepVal = NaN;
if isempty(values) || numel(values) < 2 || idx < 1 || idx > numel(values)
    return
end
if idx == 1
    return
end
prevVal = values(idx-1);
currVal = values(idx);
if ~isfinite(prevVal) || ~isfinite(currVal)
    return
end
stepVal = currVal - prevVal;
end

function window = computeParameterWindow(values)
window = [NaN, NaN];
if isempty(values)
    return
end
vals = values(isfinite(values));
if isempty(vals)
    return
end
if numel(vals) == 1
    window = [vals, vals];
else
    window = [min(vals), max(vals)];
end
end

function flag = shouldRefine(cfg)
flag = isstruct(cfg) && isfield(cfg,'refine') && isfield(cfg.refine,'onFail') && cfg.refine.onFail;
end

function tf = isTimeoutFailure(diagnostics)
tf = false;
if isempty(diagnostics) || ~isstruct(diagnostics)
    return
end
if isfield(diagnostics,'errorMessage') && ~isempty(diagnostics.errorMessage)
    msg = lower(string(diagnostics.errorMessage));
    if contains(msg, 'time limit')
        tf = true;
    end
end
end

function [sol, diagnostics] = attemptSolveWithImprovisedGuess(odeFun, bcFun, baseGuessFcn, ...
    domainMin, domainMax, stepSize, methodCfg, domainLabel)
sol = [];
lastDiag = initBranchDiagnostics(stepSize);
maxVariants = getfieldWithDefault(methodCfg,'improviseGuessAttempts',3);
for k = 1:max(1, maxVariants)
    guessFcn = improviseGuess(baseGuessFcn, domainMin, domainMax, k);
    [candidate, diag] = solveBranch(odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel);
    diag.usedContinuation = false;
    lastDiag = diag;
    if ~isempty(candidate)
        sol = candidate;
        diagnostics = diag;
        return
    end
end
diagnostics = lastDiag;
end

function guessFcnOut = improviseGuess(baseGuessFcn, domainMin, domainMax, attemptIdx)
guessFcnOut = baseGuessFcn;
if isempty(baseGuessFcn)
    return
end
seed = [];
try
    seed = baseGuessFcn(domainMin);
catch
    return
end
if isempty(seed)
    return
end
if ~isvector(seed)
    seed = seed(:,1);
end
seed = seed(:);
if ~any(isfinite(seed))
    return
end
scale = 0.2 * attemptIdx;
shift = scale * sign(seed);
shift(shift == 0) = scale;
if numel(shift) >= 3
    shift(3) = -3 * shift(3);
end
span = max(1, domainMax - domainMin);
decayRate = min(0.8, (1 + 0.3 * attemptIdx) / span);
guessFcnOut = @(x) applyGuessPerturbation(baseGuessFcn, x, domainMin, shift, decayRate);
end

function y = applyGuessPerturbation(baseGuessFcn, x, domainMin, shift, decayRate)
y = baseGuessFcn(x);
if isempty(y)
    return
end
xRow = x(:).';
decay = exp(-decayRate * (xRow - domainMin));
if isvector(y)
    y = y(:);
    y = y + shift * decay(1);
    return
end
if size(y,2) ~= numel(xRow)
    cols = min(size(y,2), numel(xRow));
    y = y(:,1:cols);
    decay = decay(1:cols);
end
y = y + shift * decay;
end

function out = enforceTimeLimit(fun, startClock, timeoutSeconds, varargin)
if toc(startClock) > timeoutSeconds
    error('displayBaseFn:timeout', 'Time limit of %.1f s exceeded for solver call.', timeoutSeconds);
end
out = fun(varargin{:});
end

function xm = midpoints(x)
if numel(x) < 2
    xm = [];
else
    xm = (x(1:end-1) + x(2:end))/2;
end
end

function [results, solutionCache, prevSolutions, branchSuccess, branchDiagnostics] = ...
    solveSweepPoint(idx, values, odeHandles, bcHandles, guesses, domainMinList, domainMaxList, gridSizeList, ...
    methodCfg, results, solutionCache, prevSolutions, odeFun, bcFun, metricEvaluators, domainLabel)

nBranches = numel(guesses);
branchSuccess = false(1, nBranches);
branchDiagnostics = repmat(initBranchDiagnostics(), 1, nBranches);

odeCurrent = odeHandles{idx};
if isempty(odeCurrent)
    odeCurrent = odeFun;
end
bcCurrent = bcHandles{idx};
if isempty(bcCurrent)
    bcCurrent = bcFun;
end

sweepValueDisplay = pickValueSafe(values, idx);

paramWindow = computeParameterWindow(values);
meshLabelDefault = getfieldWithDefault(methodCfg,'meshLabel','x');

for branchIdx = 1:nBranches
    domainMin_i = pick(domainMinList, branchIdx);
    domainMax_i = pick(domainMaxList, branchIdx);
    stepSize_i = max(3, round(pick(gridSizeList, branchIdx)));
    branchDiagnostics(branchIdx) = initBranchDiagnostics(stepSize_i);
    branchMethodCfg = methodCfg;
    branchMethodCfg.domainDisplayValue = sweepValueDisplay;
    branchMethodCfg.domainSolveRange = [domainMin_i, domainMax_i];
    branchMethodCfg.paramWindow = paramWindow;
    branchMethodCfg.paramLabel = domainLabel;
    branchMethodCfg.meshLabel = meshLabelDefault;
    if isfield(methodCfg,'groupParamName')
        branchMethodCfg.groupParamName = methodCfg.groupParamName;
    end
    if isfield(methodCfg,'groupParamValue')
        branchMethodCfg.groupParamValue = methodCfg.groupParamValue;
    end

    baseGuessFcn = guesses{branchIdx};
    nearestSol = [];
    if ~isfield(methodCfg,'disableContinuation') || ~methodCfg.disableContinuation
        nearestSol = selectNearestSolution(solutionCache(branchIdx,:), values, idx);
        if isempty(nearestSol) && ~isempty(prevSolutions{branchIdx})
            nearestSol = prevSolutions{branchIdx};
        end
    end

    [sol, usedContinuation, diag] = attemptSolveWithGuess(odeCurrent, bcCurrent, nearestSol, baseGuessFcn, ...
        domainMin_i, domainMax_i, stepSize_i, branchMethodCfg, domainLabel);
    branchDiagnostics(branchIdx) = diag;
    if isempty(sol) && usedContinuation
        [sol, usedContinuation, diag] = attemptSolveWithGuess(odeCurrent, bcCurrent, [], baseGuessFcn, ...
            domainMin_i, domainMax_i, stepSize_i, branchMethodCfg, domainLabel);
        branchDiagnostics(branchIdx) = diag;
    end
    if isempty(sol) && branchIdx == 2 && branchSuccess(1) && isTimeoutFailure(branchDiagnostics(branchIdx))
        [sol, diag] = attemptSolveWithImprovisedGuess(odeCurrent, bcCurrent, baseGuessFcn, ...
            domainMin_i, domainMax_i, stepSize_i, branchMethodCfg, domainLabel);
        branchDiagnostics(branchIdx) = diag;
    end
    branchDiagnostics(branchIdx).domainValue = branchMethodCfg.domainDisplayValue;
    if isfield(branchMethodCfg,'domainSolveRange')
        branchDiagnostics(branchIdx).domainRange = branchMethodCfg.domainSolveRange;
    end

    if isempty(sol)
        branchDiagnostics(branchIdx).stateVector = [];
        branchDiagnostics(branchIdx).solutionValue = NaN;
        branchDiagnostics(branchIdx).secondaryValue = NaN;
        printSolutionSummary([], NaN, NaN, branchDiagnostics(branchIdx), branchMethodCfg, branchIdx);
        prevSolutions{branchIdx} = [];
        solutionCache{branchIdx, idx} = [];
        results(branchIdx, idx) = blankResultStruct();
        continue
    end

    prevSolutions{branchIdx} = sol;
    solutionCache{branchIdx, idx} = sol;

    y0 = evalSolution(sol, 0);
    if ~isreal(y0), y0 = real(y0); end
    fpp0 = y0(3);
    thp0 = y0(5);

    CfVal = fpp0;
    NuVal = -thp0;
    metricContext = struct( ...
        'paramName', getfieldWithDefault(methodCfg,'paramName',''), ...
        'paramScope', getfieldWithDefault(methodCfg,'paramScope',''), ...
        'paramValue', sweepValueDisplay);
    CfVal = evaluateMetricWithFallback(metricEvaluators, 'Cf_star', sol, CfVal, metricContext);
    NuVal = evaluateMetricWithFallback(metricEvaluators, 'Nu_star', sol, NuVal, metricContext);

    duplicateTol = 1e-6;
    if usedContinuation && isDuplicateAcrossBranches(results, idx, CfVal, NuVal, branchIdx, duplicateTol)
        [solAlt, ~, diagAlt] = attemptSolveWithGuess(odeCurrent, bcCurrent, [], baseGuessFcn, ...
            domainMin_i, domainMax_i, stepSize_i, branchMethodCfg, domainLabel);
        if ~isempty(solAlt)
            y0Alt = evalSolution(solAlt, 0);
            fppAlt = y0Alt(3);
            thpAlt = y0Alt(5);
            CfAlt = fppAlt;
            NuAlt = -thpAlt;
            if ~isDuplicateAcrossBranches(results, idx, CfAlt, NuAlt, branchIdx, duplicateTol)
                sol = solAlt;
                fpp0 = fppAlt;
                thp0 = thpAlt;
                CfVal = CfAlt;
                NuVal = NuAlt;
                prevSolutions{branchIdx} = sol;
                solutionCache{branchIdx, idx} = sol;
                branchDiagnostics(branchIdx) = diagAlt;
            end
        end
    end

    results(branchIdx, idx).sol     = sol;
    results(branchIdx, idx).fpp0    = fpp0;
    results(branchIdx, idx).thp0    = thp0;
    results(branchIdx, idx).Cf_star = CfVal;
    results(branchIdx, idx).Nu_star = NuVal;

    branchDiagnostics(branchIdx).stateVector = y0(:).';
    branchDiagnostics(branchIdx).solutionValue = CfVal;
    branchDiagnostics(branchIdx).secondaryValue = NuVal;
    branchDiagnostics(branchIdx).primaryMetricLabel = 'Cf_star';
    branchDiagnostics(branchIdx).secondaryMetricLabel = 'Nu_star';

    printSolutionSummary(branchDiagnostics(branchIdx).stateVector, CfVal, NuVal, branchDiagnostics(branchIdx), branchMethodCfg, branchIdx);

    branchSuccess(branchIdx) = true;
end
end

function [sol, usedContinuation, diagnostics] = attemptSolveWithGuess(odeFun, bcFun, nearestSol, baseGuessFcn, ...
    domainMin, domainMax, stepSize, methodCfg, domainLabel)
usedContinuation = false;
guessFcn = baseGuessFcn;
if ~isempty(nearestSol)
    guessFcn = @(coord) evalSolution(nearestSol, coord);
    usedContinuation = true;
end
[sol, diagnostics] = solveBranch(odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel);
diagnostics.usedContinuation = usedContinuation;
end

function flag = isDuplicateAcrossBranches(results, idx, cfVal, nuVal, branchIdx, tol)
flag = false;
if nargin < 6 || isempty(tol)
    tol = 1e-6;
end
for other = 1:size(results,1)
    if other == branchIdx
        continue
    end
    entry = results(other, idx);
    if isempty(entry.sol)
        continue
    end
    if abs(entry.Cf_star - cfVal) < tol && abs(entry.Nu_star - nuVal) < tol
        flag = true;
        return
    end
end
end

function [values, odeHandles, bcHandles, labels, results, solutionCache, processed, prevSolutions, attemptEntries] = ...
    handleFailedSweepPoint(idx, values, odeHandles, bcHandles, labels, results, solutionCache, processed, prevSolutions, ...
    guesses, domainMinList, domainMaxList, gridSizeList, methodCfg, cfg, odeFun, bcFun, metricEvaluators, domainLabel)

attemptEntries = emptyAttemptLog();
if isempty(values)
    return
end

nBranches = size(results,1);
currIdx = min(idx, numel(values));

if isnan(values(1)) || numel(values) == 1
    [values, odeHandles, bcHandles, labels, results, solutionCache, processed] = ...
        removeValueAtIndex(values, odeHandles, bcHandles, labels, results, solutionCache, processed, currIdx);
    return
end

insertedBefore = 0;
if currIdx > 1
    [values, odeHandles, bcHandles, labels, results, solutionCache, processed, prevSolutions, insertedBefore, newAttempts] = ...
        refineBetween(values, odeHandles, bcHandles, labels, results, solutionCache, processed, prevSolutions, ...
        currIdx-1, currIdx, 'backward', cfg, guesses, domainMinList, domainMaxList, gridSizeList, ...
        methodCfg, odeFun, bcFun, nBranches, metricEvaluators, domainLabel);
    if ~isempty(newAttempts)
        attemptEntries = [attemptEntries, newAttempts]; %#ok<AGROW>
    end
    currIdx = currIdx + insertedBefore;
end

if insertedBefore == 0 && currIdx < numel(values)
    [values, odeHandles, bcHandles, labels, results, solutionCache, processed, prevSolutions, ~, newAttempts] = ...
        refineBetween(values, odeHandles, bcHandles, labels, results, solutionCache, processed, prevSolutions, ...
        currIdx, currIdx+1, 'forward', cfg, guesses, domainMinList, domainMaxList, gridSizeList, ...
        methodCfg, odeFun, bcFun, nBranches, metricEvaluators, domainLabel);
    if ~isempty(newAttempts)
        attemptEntries = [attemptEntries, newAttempts]; %#ok<AGROW>
    end
end

[values, odeHandles, bcHandles, labels, results, solutionCache, processed] = ...
    removeValueAtIndex(values, odeHandles, bcHandles, labels, results, solutionCache, processed, currIdx);
end

function [values, odeHandles, bcHandles, labels, results, solutionCache, processed, prevSolutions, insertedCount, attemptEntries] = ...
    refineBetween(values, odeHandles, bcHandles, labels, results, solutionCache, processed, prevSolutions, ...
    leftIdx, rightIdx, direction, cfg, guesses, domainMinList, domainMaxList, gridSizeList, ...
    methodCfg, odeFun, bcFun, nBranches, metricEvaluators, domainLabel)

insertedCount = 0;
attemptEntries = emptyAttemptLog();
if leftIdx < 1 || rightIdx > numel(values)
    return
end

if isfield(cfg,'refine') && isfield(cfg.refine,'numSubdiv')
    numSub = cfg.refine.numSubdiv;
else
    numSub = 5;
end

candidates = generateIntermediateValues(values(leftIdx), values(rightIdx), numSub);
if isempty(candidates)
    return
end

insertPos = rightIdx;

for k = 1:numel(candidates)
    val = candidates(k);
    [values, odeHandles, bcHandles, labels, results, solutionCache, processed, newIdx] = ...
        insertSingleValue(values, odeHandles, bcHandles, labels, results, solutionCache, processed, insertPos, val, cfg, nBranches);

    [results, solutionCache, prevSolutions, branchSuccess, branchDiagnostics] = ...
        solveSweepPoint(newIdx, values, odeHandles, bcHandles, guesses, domainMinList, domainMaxList, gridSizeList, ...
        methodCfg, results, solutionCache, prevSolutions, odeFun, bcFun, metricEvaluators, domainLabel);
    processed(newIdx) = true;
    if ~isempty(labels)
        currentLabel = labels{min(newIdx, numel(labels))};
    else
        currentLabel = '';
    end
    attemptEntries(end+1) = buildAttemptLogEntry(val, currentLabel, branchDiagnostics); %#ok<AGROW>

    if any(branchSuccess)
        insertedCount = 1;
        return
    else
        [values, odeHandles, bcHandles, labels, results, solutionCache, processed] = ...
            removeValueAtIndex(values, odeHandles, bcHandles, labels, results, solutionCache, processed, newIdx);
        if newIdx < insertPos
            insertPos = insertPos - 1;
        end
    end
end
end

function [values, odeHandles, bcHandles, labels, results, solutionCache, processed, newIdx] = ...
    insertSingleValue(values, odeHandles, bcHandles, labels, results, solutionCache, processed, insertPos, val, cfg, nBranches)

blankCol = repmat(blankResultStruct(), nBranches, 1);
values      = [values(1:insertPos-1), val, values(insertPos:end)];
odeHandles  = [odeHandles(1:insertPos-1), createHandles(cfg,'ode',val), odeHandles(insertPos:end)];
bcHandles   = [bcHandles(1:insertPos-1),  createHandles(cfg,'bc',val),  bcHandles(insertPos:end)];
labels      = [labels(1:insertPos-1),     createLabels(cfg,val),        labels(insertPos:end)];
results     = [results(:,1:insertPos-1), blankCol, results(:,insertPos:end)];
solutionCache = [solutionCache(:,1:insertPos-1), cell(nBranches,1), solutionCache(:,insertPos:end)];
processed   = [processed(1:insertPos-1), false, processed(insertPos:end)];
newIdx = insertPos;
end

function [values, odeHandles, bcHandles, labels, results, solutionCache, processed] = ...
    removeValueAtIndex(values, odeHandles, bcHandles, labels, results, solutionCache, processed, idx)

if isempty(values)
    return
end

values(idx) = [];
odeHandles(idx) = [];
bcHandles(idx) = [];
labels(idx) = [];
results(:,idx) = [];
solutionCache(:,idx) = [];
processed(idx) = [];
end

function handles = createHandles(cfg, type, vals)
handles = cell(1,numel(vals));
for k = 1:numel(vals)
    val = vals(k);
    if strcmp(type,'ode')
        if isfield(cfg,'odeFactory') && ~isempty(cfg.odeFactory)
            handles{k} = cfg.odeFactory(val);
        else
            handles{k} = [];
        end
    else
        if isfield(cfg,'bcFactory') && ~isempty(cfg.bcFactory)
            handles{k} = cfg.bcFactory(val);
        else
            handles{k} = [];
        end
    end
end
end

function labels = createLabels(cfg, vals)
labels = cell(1,numel(vals));
if isfield(cfg,'labelFcn') && ~isempty(cfg.labelFcn)
    for k = 1:numel(vals)
        labels{k} = formatLegendLabel(cfg.labelFcn(vals(k)));
    end
else
    for k = 1:numel(vals)
        labels{k} = formatNumericToken(vals(k));
    end
end
end

function s = blankResultStruct()
s = struct('sol',[],'fpp0',NaN,'thp0',NaN,'Cf_star',NaN,'Nu_star',NaN);
end

function vals = generateIntermediateValues(leftVal, rightVal, numSub)
if numSub < 1 || leftVal == rightVal
    vals = [];
    return
end
vals = linspace(leftVal, rightVal, numSub + 2);
vals = vals(2:end-1);
end

function cfg = resolveProbeConfig(sweepCfg)
cfg = struct('enabled', false, 'maxConsecutiveFails', 0);
if ~isstruct(sweepCfg)
    return
end
if isfield(sweepCfg,'probeSides') && ~isempty(sweepCfg.probeSides)
    cfg.enabled = logical(sweepCfg.probeSides);
end
if isfield(sweepCfg,'maxConsecutiveFails') && ~isempty(sweepCfg.maxConsecutiveFails)
    cfg.maxConsecutiveFails = max(0, round(sweepCfg.maxConsecutiveFails));
end
end

function state = initProbeState(cfg)
state = struct('enabled', cfg.enabled, 'mode', 'default', ...
    'leftChecked', false, 'rightChecked', false, 'rightSuccess', false, ...
    'lastProbeSide', '');
if state.enabled
    state.mode = 'probe';
end
end

function [idx, state] = selectNextSweepIndex(values, processed, cfg, state)
idx = [];
state.lastProbeSide = '';
if isempty(values) || isempty(processed)
    return
end
if ~state.enabled || ~strcmp(state.mode, 'probe')
    if strcmp(state.mode, 'right_to_left')
        remainingIdx = find(~processed);
        if ~isempty(remainingIdx)
            idx = remainingIdx(end);
        end
        return
    end
    remainingOrder = remainingSolveOrder(values, processed, cfg);
    if ~isempty(remainingOrder)
        idx = remainingOrder(1);
    end
    return
end

leftIdx = find(~processed, 1, 'first');
rightIdx = find(~processed, 1, 'last');
if isempty(leftIdx)
    return
end
if ~state.leftChecked
    idx = leftIdx;
    state.lastProbeSide = 'left';
    return
end
if ~state.rightChecked
    idx = rightIdx;
    state.lastProbeSide = 'right';
    return
end
state.mode = 'default';
remainingOrder = remainingSolveOrder(values, processed, cfg);
if ~isempty(remainingOrder)
    idx = remainingOrder(1);
end
end

function state = updateProbeStateAfterAttempt(state, values, processed)
if ~state.enabled || ~strcmp(state.mode, 'probe')
    return
end
if strcmp(state.lastProbeSide, 'left')
    state.leftChecked = true;
elseif strcmp(state.lastProbeSide, 'right')
    state.rightChecked = true;
end
if state.leftChecked && state.rightChecked
    if state.rightSuccess
        state.mode = 'right_to_left';
    else
        state.mode = 'default';
    end
end
end

function logEntries = appendSkipAttempts(logEntries, values, labels, processed, nBranches, reasonText)
if nargin < 6 || strlength(string(reasonText)) == 0
    reasonText = 'Skipped remaining sweep values.';
end
if isempty(values) || isempty(processed)
    return
end
remainingIdx = find(~processed);
if isempty(remainingIdx)
    return
end
for k = 1:numel(remainingIdx)
    idx = remainingIdx(k);
    value = pickValueSafe(values, idx);
    if isempty(labels)
        label = '';
    else
        label = labels{min(idx, numel(labels))};
    end
    diagnostics = repmat(initBranchDiagnostics(), 1, nBranches);
    for b = 1:nBranches
        diagnostics(b).status = 'fail';
        diagnostics(b).errorMessage = reasonText;
        diagnostics(b).domainValue = value;
    end
    logEntries(end+1) = buildAttemptLogEntry(value, label, diagnostics); %#ok<AGROW>
end
end

function figurePayloads = renderFigures(results, sweepValues, sweepLabels, colors, lineStyles, figConfigs, attemptLog, failRanges, successRanges, successOutliers, metricFallbacks, doPlot)
figurePayloads = {};
if isempty(figConfigs) || isempty(results)
    return
end
numSweeps = size(results,2);
if numSweeps == 0
    return
end
nBranches = size(results,1);

figurePayloads = cell(numel(figConfigs),1);
for cfgIdx = 1:numel(figConfigs)
    cfg = figConfigs(cfgIdx);
    if ~isfield(cfg,'mode') || strcmpi(cfg.mode,'profile')
        lineSets = renderProfileFigure(cfg, results, sweepLabels, colors, lineStyles, numSweeps, nBranches, sweepValues, doPlot);
        modeStr = 'profile';
    elseif strcmpi(cfg.mode,'metric')
        lineSets = renderMetricFigure(cfg, results, sweepValues, colors, lineStyles, numSweeps, nBranches, successOutliers, metricFallbacks, doPlot);
        modeStr = 'metric';
    else
        lineSets = [];
        modeStr = cfg.mode;
    end
    figurePayloads{cfgIdx} = struct( ...
        'config', cfg, ...
        'mode', modeStr, ...
        'figureId', cfg.figureId, ...
        'lineSets', {lineSets}, ...
        'attemptLog', attemptLog, ...
        'failRanges', failRanges, ...
        'successRanges', successRanges, ...
        'successOutliers', successOutliers, ...
        'branchCount', nBranches);
end
end

function lineSets = renderProfileFigure(cfg, results, sweepLabels, colors, lineStyles, numSweeps, nBranches, sweepValues, doPlot)
lineSets = struct('branchIdx',{},'sweepIdx',{},'sweepValue',{},'sweepLabel',{},'lineLabel',{},'status',{},'x',{},'y',{});
if numSweeps == 0
    return
end
domainWindow = resolveDomainWindow(cfg);
if isempty(domainWindow)
    defaultMax = 20;
    firstSol = [];
    for branchIdx = 1:nBranches
        candidate = results(branchIdx,1).sol;
        if ~isempty(candidate)
            firstSol = candidate;
            break
        end
    end
    if ~isempty(firstSol)
        defaultMax = min(firstSol.x(end), defaultMax);
    end
    domainWindow = [0, defaultMax];
end
if ~isfield(cfg,'numSamples') || isempty(cfg.numSamples)
    cfg.numSamples = 800;
end
cfg.domainWindow = domainWindow;

domainTraces = {};
codomainTraces = {};
legendEntries = {};
if doPlot
    figure(cfg.figureId); clf; hold on; grid on;
    xlabel(cfg.xlabel); ylabel(cfg.ylabel); title(cfg.title);
end

for sweepIdx = 1:numSweeps
    clr = colors(mod(sweepIdx-1, size(colors,1))+1,:);
    sweepLabel = sweepLabels{sweepIdx};
    sweepValue = pickValueSafe(sweepValues, sweepIdx);
    for branchIdx = 1:nBranches
        ls = lineStyles{mod(branchIdx-1, numel(lineStyles))+1};
        lineColor = branchColorVariant(clr, branchIdx);
        lineStruct = struct( ...
            'branchIdx', branchIdx, ...
            'sweepIdx', sweepIdx, ...
            'sweepValue', sweepValue, ...
            'sweepLabel', sweepLabel, ...
            'lineLabel', composeLegendEntry(sweepLabel, branchIdx), ...
            'status', 'no_data', ...
            'x', [], ...
            'y', []);

        sol = results(branchIdx,sweepIdx).sol;
        if isempty(sol)
            lineSets(end+1) = lineStruct; %#ok<AGROW>
            continue
        end
        domainStart = max(cfg.domainWindow(1), sol.x(1));
        domainEnd   = min(cfg.domainWindow(2), sol.x(end));
        if domainEnd <= domainStart
            lineSets(end+1) = lineStruct; %#ok<AGROW>
            continue
        end
        domainSamples = linspace(domainStart, domainEnd, cfg.numSamples);
        sampleMatrix = evalSolution(sol, domainSamples);
        values = real(sampleMatrix(cfg.componentIdx,:));
        lineStruct.status = 'success';
        lineStruct.x = domainSamples(:);
        lineStruct.y = values(:);
        lineSets(end+1) = lineStruct; %#ok<AGROW>

        if doPlot
            plot(domainSamples, values, 'LineWidth', 1.5, 'Color', lineColor, 'LineStyle', ls);
            domainTraces{end+1} = domainSamples; %#ok<AGROW>
            codomainTraces{end+1} = values; %#ok<AGROW>
            legendEntries{end+1} = lineStruct.lineLabel; %#ok<AGROW>
        end
    end
end

if doPlot
    applyDomainAxes(cfg, domainTraces, codomainTraces);
    applyNumericAxisFormatting(gca);
    if ~isempty(legendEntries)
        legendEntries = cellfun(@formatLegendLabel, legendEntries, 'UniformOutput', false);
        legend(legendEntries,'Location','best'); box on;
    end
end
end

function lineSets = renderMetricFigure(cfg, results, sweepValues, colors, lineStyles, numSweeps, nBranches, successOutliers, metricFallbacks, doPlot)
lineSets = struct('branchIdx',{},'lineLabel',{},'status',{},'x',{},'y',{},'asymptoteX',{});
if numSweeps == 0 || isempty(sweepValues)
    return
end
if numel(sweepValues) == 1 && isnan(sweepValues)
    return
end
if ~isfield(cfg,'branchIdx') || isempty(cfg.branchIdx)
    cfg.branchIdx = 1:nBranches;
end
if ~isfield(cfg,'metricField')
    cfg.metricField = 'Cf_star';
end
if ~isfield(cfg,'metricFn')
    cfg.metricFn = [];
end
if ~isfield(cfg,'xlabel') || isempty(cfg.xlabel)
    cfg.xlabel = 'sweep parameter';
end
if ~isfield(cfg,'useSharedColor')
    cfg.useSharedColor = false;
end

legendEntries = {};
yValsAll = [];
finiteSweep = sweepValues(isfinite(sweepValues));
if isempty(finiteSweep)
    return
end
metricMatrix = [];
zeroDevMask = false(1, numSweeps);
if shouldHideZeroDeviation(cfg)
    metricMatrix = buildMetricMatrix(results, cfg, numSweeps, nBranches, metricFallbacks, sweepValues);
    zeroDevMask = computeZeroDeviationMask(metricMatrix, cfg);
end

if doPlot
    figure(cfg.figureId); clf; hold on; grid on;
    xlabel(cfg.xlabel); ylabel(cfg.ylabel); title(cfg.title);
end

for idx = 1:numel(cfg.branchIdx)
    branch = cfg.branchIdx(idx);
    lineStruct = struct( ...
        'branchIdx', branch, ...
        'lineLabel', sprintf('branch %d', branch), ...
        'status', 'no_data', ...
        'x', [], ...
        'y', [], ...
        'asymptoteX', []);
    if ~isempty(metricMatrix) && idx <= size(metricMatrix, 1)
        yVals = metricMatrix(idx, :);
    else
        yVals = zeros(1, numSweeps);
        for sweepIdx = 1:numSweeps
            try
                entry = results(branch, sweepIdx);
                if ~isempty(cfg.metricFn)
                    if ~isempty(entry.sol)
                        yVals(sweepIdx) = cfg.metricFn(entry.sol);
                    else
                        yVals(sweepIdx) = NaN;
                    end
                else
                    if ~isempty(entry.sol)
                        yVals(sweepIdx) = entry.(cfg.metricField);
                    else
                        [fallbackVal, fallbackLabel] = fetchFallbackMetric(metricFallbacks, branch, pickValueSafe(sweepValues, sweepIdx), cfg.metricField);
                        if isfinite(fallbackVal)
                            yVals(sweepIdx) = fallbackVal;
                            logMetricFallbackUsage(cfg, branch, pickValueSafe(sweepValues, sweepIdx), fallbackLabel);
                        else
                            yVals(sweepIdx) = NaN;
                        end
                    end
                end
            catch
                yVals(sweepIdx) = NaN;
            end
        end
    end
    yVals = real(yVals);
    sweepValsLocal = sweepValues(:).';
    zeroDevLocal = zeroDevMask(:).';
    minLen = min([numel(yVals), numel(sweepValsLocal), numel(zeroDevLocal)]);
    if minLen < 1
        lineSets(end+1) = lineStruct; %#ok<AGROW>
        continue
    end
    yVals = yVals(1:minLen);
    sweepValsLocal = sweepValsLocal(1:minLen);
    zeroDevLocal = zeroDevLocal(1:minLen);
    validMask = isfinite(yVals) & isfinite(sweepValsLocal) & ~zeroDevLocal;
    if ~any(validMask)
        lineSets(end+1) = lineStruct; %#ok<AGROW>
        continue
    end
    lineStruct.status = 'success';
    lineStruct.x = sweepValsLocal(validMask).';
    lineStruct.y = yVals(validMask).';
    if ~shouldPlotOutliers(cfg)
        filterMode = getOutlierFilterMode(cfg);
        lineStruct = filterOutlierPoints(lineStruct, successOutliers, cfg.metricField, branch, filterMode);
    end
    % zero-deviation points already filtered via zeroDevMask
    lineStruct.asymptoteX = [];
    if isfield(cfg,'asymptote') && isstruct(cfg.asymptote) && cfg.asymptote.enabled
        [~, ~, breakIdx, asymX] = detectAsymptoticMask(lineStruct.x, lineStruct.y, getAsymptoteConfig(cfg));
        lineStruct.asymptoteX = asymX;
        if ~isempty(breakIdx)
            [lineStruct.x, lineStruct.y] = insertNaNIntoSeries(lineStruct.x, lineStruct.y, breakIdx);
        end
    end
    lineSets(end+1) = lineStruct; %#ok<AGROW>

    if doPlot
        ls = lineStyles{mod(branch-1, numel(lineStyles))+1};
        if cfg.useSharedColor
            baseColor = colors(1,:);
        else
            baseColor = colors(mod(idx-1, size(colors,1))+1,:);
        end
        branchColor = branchColorVariant(baseColor, branch);
        plotMetricLineWithAsymptote(cfg, lineStruct.x, lineStruct.y, branchColor, ls);
        legendEntries{end+1} = lineStruct.lineLabel; %#ok<AGROW>
    end
    yValsAll = [yValsAll, lineStruct.y.']; %#ok<AGROW>
end

if doPlot
    xlim([min(finiteSweep), max(finiteSweep)]);
    if ~isempty(yValsAll)
        if ~isempty(cfg.codomainWindow)
            ylim(cfg.codomainWindow);
        else
            margin = 0.05 * max(1e-6, max(yValsAll) - min(yValsAll));
            ylim([min(yValsAll)-margin, max(yValsAll)+margin]);
        end
    end
    if shouldPlotOutliers(cfg)
        plotMetricOutliers(cfg, successOutliers);
    end
    plotAsymptoteXs(collectAsymptoteX(lineSets), blendColor([0, 0, 0], [1, 1, 1], 0.5));
    applyNumericAxisFormatting(gca);
    if ~isempty(legendEntries)
        legendEntries = cellfun(@formatLegendLabel, legendEntries, 'UniformOutput', false);
        legend(legendEntries,'Location','best'); box on;
    end
end
end

function plotMetricOutliers(cfg, successOutliers)
if ~shouldPlotOutliers(cfg)
    return
end
if isempty(successOutliers) || ~isfield(cfg,'metricField')
    return
end
fieldName = cfg.metricField;
outliers = successOutliers;
if isempty(outliers)
    return
end
outlierColor = [0.85, 0.2, 0.2];
plotted = false;
for idx = 1:numel(outliers)
    entry = outliers(idx);
    if ~isfield(entry,'value')
        continue
    end
    xVal = entry.value;
    if strcmpi(fieldName, 'Nu_star') && isfield(entry,'secondaryMetric')
        yVal = entry.secondaryMetric;
    elseif isfield(entry,'primaryMetric')
        yVal = entry.primaryMetric;
    else
        continue
    end
    if ~isfinite(xVal) || ~isfinite(yVal)
        continue
    end
    h = plot(xVal, yVal, 'o', 'MarkerSize', 6, 'LineWidth', 1.2, ...
        'MarkerEdgeColor', outlierColor, 'MarkerFaceColor', 'none');
    if ~plotted
        plotted = true;
        uistack(h, 'top');
    end
end
end

function tf = shouldPlotOutliers(cfg)
tf = false;
if nargin < 1 || ~isstruct(cfg) || ~isfield(cfg,'showOutliers')
    return
end
val = cfg.showOutliers;
if islogical(val) && isscalar(val)
    tf = val;
elseif isnumeric(val) && isscalar(val)
    tf = (val ~= 0);
elseif isstring(val) && isscalar(val)
    tf = any(strcmpi(strtrim(val), ["true","1","yes","on"]));
elseif ischar(val)
    tf = any(strcmpi(strtrim(string(val)), ["true","1","yes","on"]));
end
end

function plotMetricLineWithAsymptote(cfg, xVals, yVals, baseColor, lineStyle)
asymCfg = getAsymptoteConfig(cfg);
if ~asymCfg.enabled
    plot(xVals, yVals, 'LineWidth', 1.5, 'Color', baseColor, 'LineStyle', lineStyle);
    return
end
[normalMask, asymMask, breakIdx] = detectAsymptoticMask(xVals, yVals, asymCfg);
if ~isempty(breakIdx)
    [xVals, yVals, normalMask] = insertNaNBreaks(xVals, yVals, normalMask, breakIdx);
end
plotMaskedLine(xVals, yVals, normalMask, baseColor, lineStyle);
if any(asymMask)
    asymColor = blendColor(baseColor, [1, 1, 1], asymCfg.fadeFactor);
    if ~isempty(breakIdx)
        [xAsym, yAsym, asymMask] = insertNaNBreaks(xVals, yVals, asymMask, breakIdx);
        plotMaskedLine(xAsym, yAsym, asymMask, asymColor, '--');
    else
        plotMaskedLine(xVals, yVals, asymMask, asymColor, '--');
    end
end
if asymCfg.breakLine
    plotAsymptoteBreaks(xVals, yVals, breakIdx, baseColor);
end
end

function plotMaskedLine(xVals, yVals, mask, color, style)
if isempty(xVals) || isempty(yVals)
    return
end
maskedY = yVals(:);
maskedY(~mask(:)) = NaN;
plot(xVals, maskedY, 'LineWidth', 1.5, 'Color', color, 'LineStyle', style);
end

function [xOut, yOut, maskOut] = insertNaNBreaks(xVals, yVals, mask, breakIdx)
idx = unique(breakIdx(:).');
xOut = xVals(:);
yOut = yVals(:);
maskOut = mask(:);
if isempty(idx)
    return
end
xTmp = nan(numel(xOut) + numel(idx), 1);
yTmp = nan(numel(yOut) + numel(idx), 1);
maskTmp = false(numel(maskOut) + numel(idx), 1);
srcIdx = 1;
dstIdx = 1;
for k = 1:numel(idx)
    stopIdx = idx(k);
    count = stopIdx - srcIdx + 1;
    if count > 0
        xTmp(dstIdx:dstIdx+count-1) = xOut(srcIdx:stopIdx);
        yTmp(dstIdx:dstIdx+count-1) = yOut(srcIdx:stopIdx);
        maskTmp(dstIdx:dstIdx+count-1) = maskOut(srcIdx:stopIdx);
        dstIdx = dstIdx + count;
        srcIdx = stopIdx + 1;
    end
    xTmp(dstIdx) = NaN;
    yTmp(dstIdx) = NaN;
    maskTmp(dstIdx) = false;
    dstIdx = dstIdx + 1;
end
if srcIdx <= numel(xOut)
    count = numel(xOut) - srcIdx + 1;
    xTmp(dstIdx:dstIdx+count-1) = xOut(srcIdx:end);
    yTmp(dstIdx:dstIdx+count-1) = yOut(srcIdx:end);
    maskTmp(dstIdx:dstIdx+count-1) = maskOut(srcIdx:end);
end
xOut = xTmp;
yOut = yTmp;
maskOut = maskTmp;
end

function [xOut, yOut] = insertNaNIntoSeries(xVals, yVals, breakIdx)
mask = true(size(yVals(:)));
[xOut, yOut, ~] = insertNaNBreaks(xVals, yVals, mask, breakIdx);
end

% Objective: Detect asymptotic breakpoints for metric curves.
% Purpose: Flag discontinuities based on slope flips and amplitude.
% SWOT: S-robust to noise; W-threshold tuning; O-reuse across plots; T-overfitting edge cases.
function [normalMask, asymMask, breakIdx, asymX] = detectAsymptoticMask(xVals, yVals, asymCfg)
normalMask = true(size(yVals));
asymMask = false(size(yVals));
breakIdx = [];
asymX = [];
if numel(xVals) < 3
    return
end
x = xVals(:);
y = yVals(:);
dx = diff(x);
dy = diff(y);
valid = isfinite(dx) & isfinite(dy) & (dx ~= 0);
if ~any(valid)
    return
end
slope = abs(dy(valid) ./ dx(valid));
medSlope = median(slope);
if ~isfinite(medSlope) || medSlope <= 0
    medSlope = 0;
end
threshold = asymCfg.slopeFactor * max(medSlope, eps);
edgeMask = false(numel(dx), 1);
edgeMask(valid) = abs(dy(valid) ./ dx(valid)) > threshold;
absY = abs(y);
absY = absY(isfinite(absY));
if isempty(absY)
    return
end
yLimit = prctile(absY, asymCfg.yPercentile);
asymPoints = isfinite(y) & abs(y) >= yLimit;
madVal = median(abs(absY - median(absY)));
if isfinite(madVal) && madVal > 0
    yFinite = y(isfinite(y));
    if ~isempty(yFinite)
        yMed = median(yFinite);
        asymPoints = asymPoints | (isfinite(y) & abs(y - yMed) >= (asymCfg.yMadFactor * 1.4826 * madVal));
    end
end
absDy = abs(dy);
dyValid = absDy(isfinite(absDy));
if ~isempty(dyValid)
    dyMed = median(dyValid);
    jumpMask = absDy > max(asymCfg.yJumpFactor * dyMed, eps);
    edgeMask = edgeMask | jumpMask;
end
asymMask(1:end-1) = asymMask(1:end-1) | edgeMask;
asymMask(2:end) = asymMask(2:end) | edgeMask;
asymMask = asymMask | asymPoints;
normalMask = ~asymMask;

breakIdx = [];
asymX = [];
if numel(dy) >= 2
    slopeFlip = sign(dy(1:end-1)) ~= sign(dy(2:end));
    slopeFlip = slopeFlip(:) & isfinite(dy(1:end-1)) & isfinite(dy(2:end));
    candIdx = find(slopeFlip);
    if ~isempty(candIdx)
        % break between x(k+1) and x(k+2) where slope flips
        segIdx = candIdx + 1;
        ampMask = isfinite(y(segIdx)) & (abs(y(segIdx)) >= yLimit);
        segIdx = segIdx(ampMask);
        if ~isempty(segIdx)
            [~, pickIdx] = max(abs(dy(segIdx)));
            breakIdx = segIdx(pickIdx);
            asymX = 0.5 * (x(breakIdx) + x(breakIdx + 1));
        end
    end
end
end

function cfg = getAsymptoteConfig(figCfg)
cfg = struct('enabled', false, 'slopeFactor', 10, 'yPercentile', 98, 'yMadFactor', 8, 'yJumpFactor', 10, 'fadeFactor', 0.5, 'breakLine', true);
if nargin < 1 || ~isstruct(figCfg)
    return
end
if isfield(figCfg,'asymptote') && isstruct(figCfg.asymptote)
    cfg = mergeAsymptoteConfig(cfg, figCfg.asymptote);
end
cfg.enabled = logical(cfg.enabled);
cfg.slopeFactor = max(1, cfg.slopeFactor);
cfg.yPercentile = min(100, max(0, cfg.yPercentile));
cfg.yMadFactor = max(0, cfg.yMadFactor);
cfg.yJumpFactor = max(1, cfg.yJumpFactor);
cfg.fadeFactor = min(1, max(0, cfg.fadeFactor));
cfg.breakLine = logical(cfg.breakLine);
end

function cfg = mergeAsymptoteConfig(baseCfg, updateCfg)
cfg = baseCfg;
if isfield(updateCfg,'enabled')
    cfg.enabled = updateCfg.enabled;
end
if isfield(updateCfg,'slopeFactor')
    cfg.slopeFactor = updateCfg.slopeFactor;
end
if isfield(updateCfg,'yPercentile')
    cfg.yPercentile = updateCfg.yPercentile;
end
if isfield(updateCfg,'yMadFactor')
    cfg.yMadFactor = updateCfg.yMadFactor;
end
if isfield(updateCfg,'yJumpFactor')
    cfg.yJumpFactor = updateCfg.yJumpFactor;
end
if isfield(updateCfg,'fadeFactor')
    cfg.fadeFactor = updateCfg.fadeFactor;
end
if isfield(updateCfg,'breakLine')
    cfg.breakLine = updateCfg.breakLine;
end
end

function out = blendColor(baseColor, targetColor, factor)
out = baseColor;
if numel(baseColor) ~= 3 || numel(targetColor) ~= 3
    return
end
out = baseColor + factor * (targetColor - baseColor);
out = min(1, max(0, out));
end

function plotAsymptoteBreaks(xVals, yVals, breakIdx, baseColor)
x = xVals(:);
y = yVals(:);
plotVerticalBreaks(x, y, breakIdx, baseColor);
end

function plotVerticalBreaks(xVals, yVals, breakIdx, color)
if isempty(breakIdx)
    return
end
for k = 1:numel(breakIdx)
    idx = breakIdx(k);
    if idx < 1 || idx >= numel(xVals)
        continue
    end
    if ~isfinite(xVals(idx)) || ~isfinite(xVals(idx+1))
        continue
    end
    xMid = 0.5 * (xVals(idx) + xVals(idx+1));
    yPair = yVals([idx, idx+1]);
    yPair = yPair(isfinite(yPair));
    if isempty(yPair)
        continue
    end
    yMin = min(yPair);
    yMax = max(yPair);
    plot([xMid xMid], [yMin yMax], '-.', 'LineWidth', 1.0, 'Color', color, 'HandleVisibility', 'off');
end
end

function asymX = collectAsymptoteX(lineSets)
asymX = [];
if isempty(lineSets)
    return
end
for idx = 1:numel(lineSets)
    if isfield(lineSets(idx), 'asymptoteX') && ~isempty(lineSets(idx).asymptoteX)
        asymX = [asymX; lineSets(idx).asymptoteX(:)]; %#ok<AGROW>
    end
end
end

function plotAsymptoteXs(asymX, color)
if isempty(asymX)
    return
end
ax = gca;
yl = get(ax, 'YLim');
if isempty(yl) || numel(yl) ~= 2
    return
end
uniqX = unique(asymX(isfinite(asymX)));
for k = 1:numel(uniqX)
    plot([uniqX(k) uniqX(k)], yl, '--', 'LineWidth', 1.0, 'Color', color, 'HandleVisibility', 'off');
end
end
function lineStruct = filterOutlierPoints(lineStruct, successOutliers, metricField, branchIdx, filterMode)
if isempty(successOutliers) || ~isfield(lineStruct,'x') || isempty(lineStruct.x)
    return
end
if nargin < 3 || isempty(metricField)
    metricField = '';
end
if nargin < 4
    branchIdx = NaN;
end
if nargin < 5 || isempty(filterMode)
    filterMode = 'global';
end
outlierXY = extractOutlierPoints(successOutliers, metricField, branchIdx, filterMode);
if isempty(outlierXY)
    return
end
xVals = lineStruct.x(:);
yVals = lineStruct.y(:);
keepMask = true(numel(xVals), 1);
for idx = 1:numel(xVals)
    if isOutlierPoint(xVals(idx), yVals(idx), outlierXY)
        keepMask(idx) = false;
    end
end
lineStruct.x = xVals(keepMask);
lineStruct.y = yVals(keepMask);
end

function outlierXY = extractOutlierPoints(successOutliers, metricField, branchIdx, filterMode)
outlierXY = [];
if isempty(successOutliers)
    return
end
restrictToBranch = shouldRestrictOutliers(filterMode);
for idx = 1:numel(successOutliers)
    entry = successOutliers(idx);
    if ~isfield(entry,'value')
        continue
    end
    if restrictToBranch && isfinite(branchIdx) && isfield(entry,'branchIdx') && entry.branchIdx ~= branchIdx
        continue
    end
    xVal = entry.value;
    if strcmpi(metricField, 'Nu_star') && isfield(entry,'secondaryMetric')
        yVal = entry.secondaryMetric;
    elseif isfield(entry,'primaryMetric')
        yVal = entry.primaryMetric;
    else
        continue
    end
    if isfinite(xVal) && isfinite(yVal)
        outlierXY(end+1,:) = [xVal, yVal]; %#ok<AGROW>
    end
end
end

function mode = getOutlierFilterMode(cfg)
mode = 'global';
if nargin < 1 || ~isstruct(cfg)
    return
end
if isfield(cfg,'outlierFilter') && ~isempty(cfg.outlierFilter)
    mode = char(string(cfg.outlierFilter));
end
end

function tf = shouldRestrictOutliers(filterMode)
mode = lower(strtrim(char(string(filterMode))));
tf = any(strcmp(mode, {'branch','per-branch','per_branch','perbranch'}));
end

function tf = shouldHideZeroDeviation(cfg)
tf = false;
if nargin < 1 || ~isstruct(cfg)
    return
end
if isfield(cfg,'plotZeroDeviation') && ~isempty(cfg.plotZeroDeviation)
    val = cfg.plotZeroDeviation;
    if islogical(val) && isscalar(val)
        tf = ~val;
    elseif isnumeric(val) && isscalar(val)
        tf = (val == 0);
    end
end
end

function metricMatrix = buildMetricMatrix(results, cfg, numSweeps, nBranches, metricFallbacks, sweepValues)
branchList = cfg.branchIdx;
if isempty(branchList)
    branchList = 1:nBranches;
end
metricMatrix = nan(numel(branchList), numSweeps);
for b = 1:numel(branchList)
    branch = branchList(b);
    for sweepIdx = 1:numSweeps
        try
            entry = results(branch, sweepIdx);
            if ~isempty(cfg.metricFn)
                if ~isempty(entry.sol)
                    metricMatrix(b, sweepIdx) = cfg.metricFn(entry.sol);
                end
            else
                if ~isempty(entry.sol)
                    metricMatrix(b, sweepIdx) = entry.(cfg.metricField);
                else
                    [fallbackVal, ~] = fetchFallbackMetric(metricFallbacks, branch, pickValueSafe(sweepValues, sweepIdx), cfg.metricField);
                    if isfinite(fallbackVal)
                        metricMatrix(b, sweepIdx) = fallbackVal;
                    end
                end
            end
        catch
            metricMatrix(b, sweepIdx) = NaN;
        end
    end
end
metricMatrix = real(metricMatrix);
end

function mask = computeZeroDeviationMask(metricMatrix, cfg)
mask = false(1, size(metricMatrix, 2));
if isempty(metricMatrix)
    return
end
absTol = getZeroDeviationTolerance(cfg);
for sweepIdx = 1:size(metricMatrix, 2)
    vals = metricMatrix(:, sweepIdx);
    vals = vals(isfinite(vals));
    if numel(vals) < 2
        continue
    end
    deviation = computeBranchesDeviation(vals);
    tol = max(absTol, absTol * max(1, max(abs(vals))));
    if isfinite(deviation) && deviation <= tol
        mask(sweepIdx) = true;
    end
end
end

function tol = getZeroDeviationTolerance(cfg)
tol = 1e-6;
if nargin < 1 || ~isstruct(cfg)
    return
end
if isfield(cfg,'zeroDeviationTolerance') && ~isempty(cfg.zeroDeviationTolerance)
    candidate = cfg.zeroDeviationTolerance;
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && candidate >= 0
        tol = candidate;
    end
end
end

function tf = isOutlierPoint(xVal, yVal, outlierXY)
if isempty(outlierXY)
    tf = false;
    return
end
tol = 1e-9 * max(1, max(abs([xVal, yVal])));
dx = abs(outlierXY(:,1) - xVal);
dy = abs(outlierXY(:,2) - yVal);
tf = any(dx <= tol & dy <= tol);
end

function applyDomainAxes(cfg, domainTraces, codomainTraces)
codomainMin = inf; codomainMax = -inf;
for k = 1:numel(codomainTraces)
    vals = codomainTraces{k};
    if isempty(vals)
        continue
    end
    codomainMin = min(codomainMin, min(vals));
    codomainMax = max(codomainMax, max(vals));
end
if ~isfinite(codomainMin) || ~isfinite(codomainMax)
    codomainMin = 0; codomainMax = 1;
end

window = resolveDomainWindow(cfg);
if isempty(window)
    window = [0, 1];
end
xlim(window);
codomainWindow = resolveCodomainWindow(cfg);
if ~isempty(codomainWindow)
    ylim(codomainWindow);
else
    margin = 0.05 * max(1e-6, codomainMax - codomainMin);
    ylim([codomainMin - margin, codomainMax + margin]);
end
end

function applyNumericAxisFormatting(axHandle)
if nargin < 1 || isempty(axHandle)
    axHandle = gca;
end
try
    xTicks = get(axHandle, 'XTick');
    yTicks = get(axHandle, 'YTick');
    if ~isempty(xTicks)
        xLabels = arrayfun(@formatNumericToken, xTicks, 'UniformOutput', false);
        set(axHandle, 'XTickLabel', xLabels);
    end
    if ~isempty(yTicks)
        yLabels = arrayfun(@formatNumericToken, yTicks, 'UniformOutput', false);
        set(axHandle, 'YTickLabel', yLabels);
    end
catch
end
end

function entry = composeLegendEntry(sweepLabel, branchIdx)
if isempty(sweepLabel)
    entry = sprintf('branch %d', branchIdx);
else
    entry = sprintf('%s, branch %d', sweepLabel, branchIdx);
end
end

function out = formatLegendLabel(label)
if isstring(label)
    txt = char(label);
elseif ischar(label)
    txt = label;
else
    txt = char(string(label));
end
if isempty(txt)
    out = txt;
    return
end
out = replaceNumericTokens(txt);
end

function out = replaceNumericTokens(txt)
out = txt;
[starts, ends, matches] = regexp(txt, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'start', 'end', 'match');
if isempty(matches)
    return
end
for idx = numel(matches):-1:1
    numVal = str2double(matches{idx});
    if ~isfinite(numVal)
        continue
    end
    formatted = formatNumericToken(numVal);
    out = [out(1:starts(idx)-1), formatted, out(ends(idx)+1:end)];
end
end

function entry = buildAttemptLogEntry(value, label, diagnostics)
nBranches = numel(diagnostics);
entry = struct();
entry.value = value;
entry.label = label;
entry.branchStatus = repmat({''}, 1, nBranches);
entry.stepSizes = nan(1, nBranches);
entry.avgMeshSteps = nan(1, nBranches);
entry.maxResiduals = nan(1, nBranches);
entry.absoluteErrors = nan(1, nBranches);
entry.relativeErrors = nan(1, nBranches);
entry.errorMessages = repmat({''}, 1, nBranches);
entry.meshPoints = nan(1, nBranches);
entry.usedContinuation = false(1, nBranches);
entry.solutionValues = nan(1, nBranches);
entry.secondaryValues = nan(1, nBranches);
entry.primaryRates = nan(1, nBranches);
entry.secondaryRates = nan(1, nBranches);
entry.stateVectors = repmat({[]}, 1, nBranches);
entry.sweepSteps = nan(1, nBranches);
entry.initialGuessErrors = nan(1, nBranches);
entry.iterations = nan(1, nBranches);
entry.consoleLogs = repmat({''}, 1, nBranches);
entry.primaryMetricLabels = repmat({''}, 1, nBranches);
entry.secondaryMetricLabels = repmat({''}, 1, nBranches);
entry.initialGuessMesh = repmat({[]}, 1, nBranches);
entry.initialGuessProfiles = repmat({[]}, 1, nBranches);
entry.allBranchesSuccess = false;
entry.branchesDeviation = NaN;
entry.branchesDeviationPrimary = NaN;
entry.branchesDeviationSecondary = NaN;

for idx = 1:nBranches
    diag = diagnostics(idx);
    if isfield(diag,'status') && ~isempty(diag.status)
        entry.branchStatus{idx} = diag.status;
    end
    if isfield(diag,'requestedStep') && ~isempty(diag.requestedStep)
        entry.stepSizes(idx) = diag.requestedStep;
    end
    if isfield(diag,'avgMeshStep') && ~isempty(diag.avgMeshStep)
        entry.avgMeshSteps(idx) = diag.avgMeshStep;
    end
    if isfield(diag,'maxResidual') && ~isempty(diag.maxResidual)
        entry.maxResiduals(idx) = diag.maxResidual;
    end
    if isfield(diag,'maxResidual') && ~isempty(diag.maxResidual) && isfinite(diag.maxResidual)
        entry.absoluteErrors(idx) = diag.maxResidual;
        if isfield(diag,'solutionValue') && ~isempty(diag.solutionValue) && isfinite(diag.solutionValue)
            denom = max(eps, abs(diag.solutionValue));
            entry.relativeErrors(idx) = diag.maxResidual / denom;
        end
    end
    if isfield(diag,'errorMessage') && ~isempty(diag.errorMessage)
        entry.errorMessages{idx} = diag.errorMessage;
    end
    if isfield(diag,'meshPoints') && ~isempty(diag.meshPoints)
        entry.meshPoints(idx) = diag.meshPoints;
    end
    if isfield(diag,'usedContinuation') && ~isempty(diag.usedContinuation)
        entry.usedContinuation(idx) = logical(diag.usedContinuation);
    end
    if isfield(diag,'solutionValue') && ~isempty(diag.solutionValue)
        entry.solutionValues(idx) = diag.solutionValue;
    end
    if isfield(diag,'secondaryValue') && ~isempty(diag.secondaryValue)
        entry.secondaryValues(idx) = diag.secondaryValue;
    end
    if isfield(diag,'primaryMetricLabel') && ~isempty(diag.primaryMetricLabel)
        entry.primaryMetricLabels{idx} = diag.primaryMetricLabel;
    end
    if isfield(diag,'secondaryMetricLabel') && ~isempty(diag.secondaryMetricLabel)
        entry.secondaryMetricLabels{idx} = diag.secondaryMetricLabel;
    end
    if isfield(diag,'sweepStep') && ~isempty(diag.sweepStep)
        entry.sweepSteps(idx) = diag.sweepStep;
    end
    if isfield(diag,'initialGuessError') && ~isempty(diag.initialGuessError)
        entry.initialGuessErrors(idx) = diag.initialGuessError;
    end
    if isfield(diag,'iterations') && ~isempty(diag.iterations)
        entry.iterations(idx) = diag.iterations;
    end
    if isfield(diag,'consoleLog') && ~isempty(diag.consoleLog)
        entry.consoleLogs{idx} = diag.consoleLog;
    end
    if isfield(diag,'stateVector')
        entry.stateVectors{idx} = diag.stateVector;
    end
    if isfield(diag,'initialGuessMesh') && ~isempty(diag.initialGuessMesh)
        entry.initialGuessMesh{idx} = diag.initialGuessMesh;
    end
    if isfield(diag,'initialGuessProfile') && ~isempty(diag.initialGuessProfile)
        entry.initialGuessProfiles{idx} = diag.initialGuessProfile;
    end
end
entry.allBranchesSuccess = all(strcmpi(entry.branchStatus, 'success'));
entry.branchesDeviation = computeBranchesDeviation(entry.solutionValues);
entry.branchesDeviationPrimary = computeBranchesDeviation(entry.solutionValues);
entry.branchesDeviationSecondary = computeBranchesDeviation(entry.secondaryValues);
end

function logStruct = emptyAttemptLog()
logStruct = struct('value',{},'label',{},'branchStatus',{}, ...
    'stepSizes',{},'avgMeshSteps',{},'maxResiduals',{}, ...
    'absoluteErrors',{},'relativeErrors',{},'errorMessages',{},'meshPoints',{},'usedContinuation',{}, ...
    'solutionValues',{},'secondaryValues',{},'primaryRates',{},'secondaryRates',{}, ...
    'primaryMetricLabels',{},'secondaryMetricLabels',{}, ...
    'sweepSteps',{},'initialGuessErrors',{},'iterations',{}, ...
    'consoleLogs',{},'stateVectors',{},'initialGuessMesh',{},'initialGuessProfiles',{}, ...
    'allBranchesSuccess',{},'branchesDeviation',{},'branchesDeviationPrimary',{},'branchesDeviationSecondary',{});
end

function attemptLog = computeAttemptRates(attemptLog)
if isempty(attemptLog) || numel(attemptLog) < 2
    return
end
for idx = 2:numel(attemptLog)
    prev = attemptLog(idx-1);
    curr = attemptLog(idx);
    deltaSweep = curr.value - prev.value;
    if ~isfinite(deltaSweep) || deltaSweep == 0
        attemptLog(idx) = curr;
        continue
    end
    nBranches = max([numel(curr.solutionValues), numel(prev.solutionValues), numel(curr.secondaryValues), numel(prev.secondaryValues)]);
    primaryRates = nan(1, nBranches);
    secondaryRates = nan(1, nBranches);
    for b = 1:nBranches
        if isBranchSuccess(curr, b) && isBranchSuccess(prev, b)
            currPrimary = pickAttemptValue(curr.solutionValues, b);
            prevPrimary = pickAttemptValue(prev.solutionValues, b);
            if isfinite(currPrimary) && isfinite(prevPrimary)
                primaryRates(b) = (currPrimary - prevPrimary) / deltaSweep;
            end
            currSecondary = pickAttemptValue(curr.secondaryValues, b);
            prevSecondary = pickAttemptValue(prev.secondaryValues, b);
            if isfinite(currSecondary) && isfinite(prevSecondary)
                secondaryRates(b) = (currSecondary - prevSecondary) / deltaSweep;
            end
        end
    end
    curr.primaryRates = primaryRates;
    curr.secondaryRates = secondaryRates;
    attemptLog(idx) = curr;
end
end

function [values, sweepLabels, results, attemptLog] = applySmartRateFilter(values, sweepLabels, results, attemptLog, sweepCfg)
if isempty(values) || isempty(attemptLog) || ~isfield(sweepCfg,'smartRateEnabled') || ~sweepCfg.smartRateEnabled
    return
end
if ~isfield(sweepCfg,'smartRateConfig') || isempty(sweepCfg.smartRateConfig) || ~isstruct(sweepCfg.smartRateConfig)
    return
end
cfg = sweepCfg.smartRateConfig;
if ~isfield(cfg,'min') || ~isfield(cfg,'max') || ~isfield(cfg,'count') ...
        || ~isfield(cfg,'minRate') || ~isfield(cfg,'maxRate')
    return
end
minVal = cfg.min;
maxVal = cfg.max;
targetCount = round(cfg.count);
minRate = cfg.minRate;
maxRate = cfg.maxRate;
if targetCount < 1 || maxRate < 0 || minRate < 0
    return
end
if maxVal < minVal
    tmp = minVal;
    minVal = maxVal;
    maxVal = tmp;
end

tol = 1e-7;
[candidateIdx, candidateValues] = collectRateCandidates(attemptLog, minVal, maxVal);
if isempty(candidateIdx)
    warning('displayBaseFn:SmartRateNoMatches', ...
        'No sweep points met smart rate bounds; keeping all points.');
    return
end
if targetCount > numel(candidateIdx)
    warning('displayBaseFn:SmartRateTooFew', ...
        'Only %d points eligible for rate filtering (target %d).', numel(candidateIdx), targetCount);
    targetCount = numel(candidateIdx);
end
pathIdx = findRatePath(attemptLog, candidateIdx, targetCount, minRate, maxRate);
if isempty(pathIdx)
    fallbackCount = targetCount;
    while isempty(pathIdx) && fallbackCount > 1
        fallbackCount = fallbackCount - 1;
        pathIdx = findRatePath(attemptLog, candidateIdx, fallbackCount, minRate, maxRate);
    end
    if isempty(pathIdx)
        warning('displayBaseFn:SmartRateNoPath', ...
            'No feasible rate path found; keeping the first eligible point only.');
        pathIdx = 1;
    else
        warning('displayBaseFn:SmartRateReduced', ...
            'Reduced target count to %d to satisfy rate bounds.', fallbackCount);
    end
end
selectedVals = candidateValues(pathIdx);
selectedIdx = [];
for k = 1:numel(selectedVals)
    matchIdx = find(abs(values - selectedVals(k)) <= tol, 1, 'first');
    if ~isempty(matchIdx)
        selectedIdx(end+1) = matchIdx; %#ok<AGROW>
    end
end
selectedIdx = unique(selectedIdx, 'stable');
if isempty(selectedIdx)
    warning('displayBaseFn:SmartRateNoMatches', ...
        'No sweep points met smart rate bounds; keeping all points.');
    return
end

values = values(selectedIdx);
if ~isempty(sweepLabels)
    sweepLabels = sweepLabels(selectedIdx);
end
if ~isempty(results) && size(results,2) >= max(selectedIdx)
    results = results(:, selectedIdx);
end
attemptLog = filterAttemptLogByValues(attemptLog, values, tol);
attemptLog = computeAttemptRates(attemptLog);
end

function [idxList, values] = collectRateCandidates(attemptLog, minVal, maxVal)
idxList = [];
values = [];
for idx = 1:numel(attemptLog)
    entry = attemptLog(idx);
    if ~isfield(entry,'allBranchesSuccess') || ~entry.allBranchesSuccess
        continue
    end
    if ~isfinite(entry.value) || entry.value < minVal || entry.value > maxVal
        continue
    end
    if isempty(entry.solutionValues) || isempty(entry.secondaryValues)
        continue
    end
    idxList(end+1,1) = idx; %#ok<AGROW>
    values(end+1,1) = entry.value; %#ok<AGROW>
end
if isempty(idxList)
    return
end
[values, order] = sort(values);
idxList = idxList(order);
end

function pathIdx = findRatePath(attemptLog, idxList, targetCount, minRate, maxRate)
n = numel(idxList);
if targetCount <= 1
    pathIdx = 1;
    return
end
dp = -ones(targetCount, n);
for j = 1:n
    dp(1, j) = 0;
end
for k = 2:targetCount
    for j = 2:n
        for i = 1:j-1
            if dp(k-1, i) < 0
                continue
            end
            if rateEdgeOk(attemptLog(idxList(i)), attemptLog(idxList(j)), minRate, maxRate)
                dp(k, j) = i;
                break
            end
        end
    end
end
endIdx = find(dp(targetCount, :) >= 0, 1, 'first');
if isempty(endIdx)
    pathIdx = [];
    return
end
pathIdx = zeros(targetCount, 1);
pathIdx(targetCount) = endIdx;
for k = targetCount:-1:2
    pathIdx(k-1) = dp(k, pathIdx(k));
end
end

function ok = rateEdgeOk(prevEntry, currEntry, minRate, maxRate)
deltaSweep = currEntry.value - prevEntry.value;
if ~isfinite(deltaSweep) || deltaSweep == 0
    ok = false;
    return
end
nBranches = max([numel(prevEntry.solutionValues), numel(currEntry.solutionValues), ...
    numel(prevEntry.secondaryValues), numel(currEntry.secondaryValues)]);
for b = 1:nBranches
    prevPrimary = pickAttemptValue(prevEntry.solutionValues, b);
    currPrimary = pickAttemptValue(currEntry.solutionValues, b);
    prevSecondary = pickAttemptValue(prevEntry.secondaryValues, b);
    currSecondary = pickAttemptValue(currEntry.secondaryValues, b);
    if ~isfinite(prevPrimary) || ~isfinite(currPrimary) || ~isfinite(prevSecondary) || ~isfinite(currSecondary)
        ok = false;
        return
    end
    primaryRate = abs((currPrimary - prevPrimary) / deltaSweep);
    secondaryRate = abs((currSecondary - prevSecondary) / deltaSweep);
    if primaryRate < minRate || primaryRate > maxRate || secondaryRate < minRate || secondaryRate > maxRate
        ok = false;
        return
    end
end
ok = true;
end

function filtered = filterAttemptLogByValues(attemptLog, values, tol)
filtered = emptyAttemptLog();
if isempty(attemptLog) || isempty(values)
    return
end
for idx = 1:numel(attemptLog)
    entry = attemptLog(idx);
    if ~isfinite(entry.value)
        continue
    end
    if any(abs(values - entry.value) <= tol)
        filtered(end+1) = entry; %#ok<AGROW>
    end
end
end

function tf = isBranchSuccess(entry, branchIdx)
tf = false;
if ~isfield(entry,'branchStatus') || branchIdx > numel(entry.branchStatus)
    return
end
status = entry.branchStatus{branchIdx};
if isempty(status)
    return
end
tf = strcmpi(status, 'success');
end

function value = pickAttemptValue(values, idx)
if isempty(values)
    value = NaN;
    return
end
if numel(values) == 1
    value = values;
    return
end
value = values(min(idx, numel(values)));
end
function deviation = computeBranchesDeviation(values)
deviation = NaN;
if isempty(values)
    return
end
vals = values(isfinite(values));
if numel(vals) < 2
    return
end
numVals = numel(vals);
pairs = [];
for i = 1:numVals-1
    for j = i+1:numVals
        pairs(end+1) = abs(vals(i) - vals(j)); %#ok<AGROW>
    end
end
if isempty(pairs)
    return
end
deviation = mean(pairs);
end

function diagnostics = assignSweepStep(diagnostics, sweepStep)
for idx = 1:numel(diagnostics)
    diagnostics(idx).sweepStep = sweepStep;
end
end

function diag = initBranchDiagnostics(stepSize)
if nargin < 1
    stepSize = NaN;
end
diag = struct( ...
    'status', 'not_run', ...
    'requestedStep', stepSize, ...
    'avgMeshStep', NaN, ...
    'maxResidual', NaN, ...
    'errorMessage', '', ...
    'meshPoints', NaN, ...
    'usedContinuation', false, ...
    'warningId', '', ...
    'warningMsg', '', ...
    'attempts', 0, ...
    'solutionValue', NaN, ...
    'secondaryValue', NaN, ...
    'primaryMetricLabel', '', ...
    'secondaryMetricLabel', '', ...
    'sweepStep', NaN, ...
    'initialGuessError', NaN, ...
    'iterations', NaN, ...
    'consoleLog', '', ...
    'stateVector', [], ...
    'domainValue', NaN, ...
    'domainRange', [NaN, NaN], ...
    'initialGuessMesh', [], ...
    'initialGuessProfile', []);
end

function cfg = normalizeMethodSpec(cfg, domainMinList, domainMaxList)
if nargin < 1 || isempty(cfg)
    cfg = struct();
end
if nargin < 2
    domainMinList = 0;
end
if nargin < 3
    domainMaxList = 1;
end
if ~isfield(cfg,'solver') || strlength(string(cfg.solver)) == 0
    cfg.solver = 'bvp4c';
end
span = abs(pick(domainMaxList,1) - pick(domainMinList,1));
if ~isfinite(span) || span <= 0
    span = 1;
end
cfg.displayName = getfieldWithDefault(cfg,'displayName',upper(string(cfg.solver)));
cfg.relativeTolerance = getfieldWithDefault(cfg,'relativeTolerance',1e-6);
cfg.absoluteTolerance = getfieldWithDefault(cfg,'absoluteTolerance',1e-8);
cfg.maxMeshPoints = getfieldWithDefault(cfg,'maxMeshPoints',60000);
cfg.maxSolverAttempts = getfieldWithDefault(cfg,'maxSolverAttempts',10);
cfg.timeLimitSeconds = getfieldWithDefault(cfg,'timeLimitSeconds',10);
cfg.attemptTimeLimit = getfieldWithDefault(cfg,'attemptTimeLimit',min(10, cfg.timeLimitSeconds));
cfg.maxStepCount = getfieldWithDefault(cfg,'maxStepCount',20000);
cfg.initialStepSize = getfieldWithDefault(cfg,'initialStepSize',span/200);
cfg.minStepSize = getfieldWithDefault(cfg,'minStepSize',cfg.initialStepSize/20);
cfg.maxStepSize = getfieldWithDefault(cfg,'maxStepSize',span/10);
cfg.stepSafetyFactor = getfieldWithDefault(cfg,'stepSafetyFactor',0.9);
cfg.stepGrowthLimit = getfieldWithDefault(cfg,'stepGrowthLimit',4);
cfg.stepShrinkLimit = getfieldWithDefault(cfg,'stepShrinkLimit',0.25);
cfg.newtonTolerance = getfieldWithDefault(cfg,'newtonTolerance',1e-8);
cfg.newtonMaxIterations = getfieldWithDefault(cfg,'newtonMaxIterations',12);
cfg.jacobianPerturbation = getfieldWithDefault(cfg,'jacobianPerturbation',1e-6);
cfg.shootingTolerance = getfieldWithDefault(cfg,'shootingTolerance',1e-6);
cfg.shootingMaxIterations = getfieldWithDefault(cfg,'shootingMaxIterations',10);
cfg.shootingPerturbation = getfieldWithDefault(cfg,'shootingPerturbation',1e-6);
cfg.meshLabel = getfieldWithDefault(cfg,'meshLabel','x');
if ~isfield(cfg,'freeInitialIndices') || isempty(cfg.freeInitialIndices)
    cfg.freeInitialIndices = [];
else
    cfg.freeInitialIndices = cfg.freeInitialIndices(:).';
end
if ~isfield(cfg,'initialStateProjector')
    cfg.initialStateProjector = [];
end
if ~isfield(cfg,'activeResidualIndices') || isempty(cfg.activeResidualIndices)
    cfg.activeResidualIndices = [];
else
    cfg.activeResidualIndices = cfg.activeResidualIndices(:).';
end
end

function val = getfieldWithDefault(s, fieldName, defaultVal)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    val = s.(fieldName);
else
    val = defaultVal;
end
end

function val = pickValueSafe(vals, idx)
if isempty(vals)
    val = NaN;
    return
end
if numel(vals) == 1
    val = vals;
    return
end
val = vals(min(idx, numel(vals)));
end

function value = evaluateMetricWithFallback(metricEvaluators, fieldName, sol, fallback, context)
value = fallback;
if nargin < 1 || isempty(metricEvaluators) || ~isstruct(metricEvaluators)
    return
end
if ~isfield(metricEvaluators, fieldName)
    return
end
fn = metricEvaluators.(fieldName);
if isempty(fn) || ~isa(fn,'function_handle')
    return
end
try
    if nargin >= 5 && ~isempty(context)
        nArgs = nargin(fn);
        if nArgs < 0 || nArgs >= 2
            candidate = fn(sol, context);
        else
            candidate = fn(sol);
        end
    else
        candidate = fn(sol);
    end
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        value = candidate;
    end
catch
    value = fallback;
end
end

function printSolutionSummary(stateVector, cfValue, nuValue, diagInfo, methodMeta, branchIdx)
statusTxt = 'success';
reasonTxt = '';
attemptValueTxt = '';
if nargin < 6
    branchIdx = NaN;
end
if nargin >= 4 && isstruct(diagInfo)
    if isfield(diagInfo,'status') && ~isempty(diagInfo.status)
        statusTxt = diagInfo.status;
    end
    if isfield(diagInfo,'errorMessage') && ~isempty(diagInfo.errorMessage)
        reasonTxt = diagInfo.errorMessage;
    end
end
attemptTokens = buildAttemptTokens(methodMeta, diagInfo);
if ~isempty(attemptTokens)
    attemptValueTxt = sprintf('Attempted %s;', strjoin(attemptTokens, '; '));
end
if isempty(reasonTxt)
    reasonTxt = 'n/a';
end
logIdx = fetchNextLogIndex();
emitRunLog('(%d) %s;\n', logIdx, datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'));
if isfinite(branchIdx)
    emitRunLog('Branch: %d;\n', branchIdx);
end
emitRunLog('Status: %s;\n', statusTxt);
if ~isempty(attemptValueTxt)
    emitRunLog('%s\n', attemptValueTxt);
end
emitRunLog('Reason: %s;\n', reasonTxt);

if isfield(diagInfo,'initialGuessMesh') && ~isempty(diagInfo.initialGuessMesh)
    emitRunLog('Initial guess mesh: %s\n', formatGuessArrayForExport(diagInfo.initialGuessMesh));
end
if isfield(diagInfo,'initialGuessProfile') && ~isempty(diagInfo.initialGuessProfile)
    emitRunLog('Initial guess profile: %s\n', formatGuessArrayForExport(diagInfo.initialGuessProfile));
end
if isempty(stateVector)
    emitRunLog('Result: N/A\n');
else
    vec = stateVector(:).';
    formattedVals = strjoin(arrayfun(@formatNumericToken, vec, 'UniformOutput', false), ', ');
    emitRunLog('Result: [%s]\n', formattedVals);
end
if nargin < 2 || ~isfinite(cfValue)
    emitRunLog('localSkinFriction: N/A\n');
else
    emitRunLog('localSkinFriction: %s\n', formatNumericToken(cfValue));
end
if nargin < 3 || ~isfinite(nuValue)
    emitRunLog('nusseltNumber: N/A\n');
else
    emitRunLog('nusseltNumber: %s\n', formatNumericToken(nuValue));
end
emitRunLog('=======================\n');
end

function logText = formatSolverConsoleLog(rawLog, attemptNum, domainMin, domainMax, domainLabel, methodMeta)
logBody = strtrim(rawLog);
if isempty(logBody)
    logText = '';
    return
end
if nargin < 5 || isempty(domainLabel)
    domainLabel = 'parameter';
end
paramLabelClean = sanitizeAxisLabelLocal(domainLabel);
meshLabel = 'x';
paramWindowTxt = '';
meshWindowTxt = '';
meshSolvedTxt = '';
groupTxt = '';

if nargin >= 6 && ~isempty(methodMeta) && isstruct(methodMeta)
    if isfield(methodMeta,'paramLabel') && ~isempty(methodMeta.paramLabel)
        paramLabelClean = sanitizeAxisLabelLocal(methodMeta.paramLabel);
    end
    if isfield(methodMeta,'paramWindow')
        window = methodMeta.paramWindow;
        if numel(window) >= 2 && all(isfinite(window))
            valueSuffix = '';
            if isfield(methodMeta,'domainDisplayValue') && isfinite(methodMeta.domainDisplayValue)
                valueSuffix = sprintf(' (%s=%s)', paramLabelClean, formatNumericToken(methodMeta.domainDisplayValue));
            end
            paramWindowTxt = sprintf('%s window: [%s, %s]%s;', paramLabelClean, formatNumericToken(window(1)), formatNumericToken(window(2)), valueSuffix);
        elseif isfield(methodMeta,'domainDisplayValue') && isfinite(methodMeta.domainDisplayValue)
            paramWindowTxt = sprintf('%s=%s;', paramLabelClean, formatNumericToken(methodMeta.domainDisplayValue));
        end
    elseif isfield(methodMeta,'domainDisplayValue') && isfinite(methodMeta.domainDisplayValue)
        paramWindowTxt = sprintf('%s=%s;', paramLabelClean, formatNumericToken(methodMeta.domainDisplayValue));
    end
    if isfield(methodMeta,'meshLabel') && ~isempty(methodMeta.meshLabel)
        meshLabel = sanitizeAxisLabelLocal(methodMeta.meshLabel);
    end
    if isfield(methodMeta,'solutionRange') && numel(methodMeta.solutionRange) >= 2
        solveBounds = methodMeta.solutionRange;
        if all(isfinite(solveBounds))
            meshSolvedTxt = sprintf('%s mesh solved range [%s, %s];', meshLabel, formatNumericToken(solveBounds(1)), formatNumericToken(solveBounds(end)));
        end
    end
    if isfield(methodMeta,'groupParamName') && ~isempty(methodMeta.groupParamName)
        groupLabel = sanitizeAxisLabelLocal(methodMeta.groupParamName);
        if isfield(methodMeta,'groupParamValue') && isfinite(methodMeta.groupParamValue)
            groupTxt = sprintf('%s=%s;', groupLabel, formatNumericToken(methodMeta.groupParamValue));
        else
            groupTxt = sprintf('%s group;', groupLabel);
        end
    end
end

meshWindowTxt = sprintf('%s mesh window: [%s, %s];', meshLabel, formatNumericToken(domainMin), formatNumericToken(domainMax));

headerLines = {
    '======================='
    sprintf('%s;', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'))
    sprintf('Attempt: %d;', attemptNum)
    };
if ~isempty(groupTxt), headerLines{end+1} = groupTxt; end %#ok<AGROW>
if ~isempty(paramWindowTxt), headerLines{end+1} = paramWindowTxt; end %#ok<AGROW>
headerLines{end+1} = meshWindowTxt;
if ~isempty(meshSolvedTxt), headerLines{end+1} = meshSolvedTxt; end %#ok<AGROW>

headerText = strjoin(headerLines, newline);
logText = sprintf('%s\n%s\n=======================', headerText, logBody);
end

function err = computeInitialGuessError(solStruct, guessFcn)
err = NaN;
if nargin < 2 || isempty(solStruct) || isempty(guessFcn)
    return
end
try
    guessVals = evaluateGuessOnMesh(guessFcn, solStruct.x);
    if isempty(guessVals)
        return
    end
    solVals = solStruct.y;
    minStates = min(size(solVals,1), size(guessVals,1));
    if minStates == 0
        return
    end
    diffMat = solVals(1:minStates,:) - guessVals(1:minStates,:);
    err = sqrt(mean(diffMat(:).^2));
catch
    err = NaN;
end
end

function txt = formatGuessArrayForExport(value)
if isempty(value)
    txt = '';
    return
end
if isnumeric(value)
    txt = formatNumericArray(value);
    return
end
if iscell(value)
    txt = formatGuessArrayForExport(value{1});
    return
end
txt = char(string(value));
end

function txt = buildMethodConsoleLog(meshPoints, residualValue, odeEvalCount, bcEvalCount)
lines = {
    sprintf('The solution was obtained on a mesh of %d points.', meshPoints);
    sprintf('The maximum residual is  %s.', formatNumericToken(residualValue));
    sprintf('There were %d calls to the ODE function.', odeEvalCount);
    sprintf('There were %d calls to the BC function.', bcEvalCount);
    };
txt = strjoin(lines, newline);
end

function window = resolveDomainWindow(cfg)
window = [];
if isfield(cfg,'domainWindow') && ~isempty(cfg.domainWindow)
    window = cfg.domainWindow;
end
end

function window = resolveCodomainWindow(cfg)
window = [];
if isfield(cfg,'codomainWindow') && ~isempty(cfg.codomainWindow)
    window = cfg.codomainWindow;
end
end

function label = sanitizeAxisLabelLocal(value)
txt = char(string(value));
txt = strrep(txt, '\', '');
txt = regexprep(txt, '[\{\}^_]', '');
txt = strtrim(txt);
if isempty(txt)
    label = 'x';
else
    label = txt;
end
end

function guessVals = evaluateGuessOnMesh(guessFcn, mesh)
guessVals = [];
if isempty(guessFcn) || isempty(mesh)
    return
end
numPts = numel(mesh);
sample = guessFcn(mesh(1));
if isempty(sample)
    guessVals = [];
    return
end
stateCount = numel(sample);
guessVals = zeros(stateCount, numPts);
guessVals(:,1) = sample(:);
for idx = 2:numPts
    col = guessFcn(mesh(idx));
    colVec = col(:);
    if numel(colVec) < stateCount
        colVec(end+1:stateCount,1) = colVec(end); %#ok<AGROW>
    elseif numel(colVec) > stateCount
        colVec = colVec(1:stateCount);
    end
    guessVals(:,idx) = colVec;
end
end

function iter = extractIterationCount(solStruct)
iter = NaN;
if ~isstruct(solStruct) || ~isfield(solStruct,'stats') || isempty(solStruct.stats)
    return
end
stats = solStruct.stats;
candidateFields = {'niter','nIter','iterations','Iterations','nODEevals','nOdeEvals','nfevals','nFevals'};
for k = 1:numel(candidateFields)
    name = candidateFields{k};
    if isfield(stats, name) && ~isempty(stats.(name))
        iter = stats.(name);
        return
    end
end
end

function avgStep = averageMeshSpacing(mesh)
if nargin == 0 || numel(mesh) < 2
    avgStep = NaN;
else
    avgStep = (mesh(end) - mesh(1)) / (numel(mesh) - 1);
end
end

function maxRes = extractMaxResidual(sol)
maxRes = NaN;
if ~isstruct(sol) || ~isfield(sol,'stats') || isempty(sol.stats)
    return
end
stats = sol.stats;
if ~isstruct(stats)
    return
end
candidateFields = {'maxres','maxRes','maxresidual','maxResidual','maxresiduals','maxResid'};
for k = 1:numel(candidateFields)
    name = candidateFields{k};
    if isfield(stats, name) && ~isempty(stats.(name))
        maxRes = stats.(name);
        return
    end
end
end

function yq = evalSolution(sol, xq)
% evalSolution  Evaluate BVP solution at requested points with a deval fallback.
persistent canUseDeval
if isempty(canUseDeval)
    canUseDeval = true;
end

if canUseDeval && isstruct(sol) && isfield(sol,'solver')
    try
        yq = deval(sol, xq);
        return
    catch err
        msg = lower(err.message);
        if contains(msg, 'input arguments') || contains(msg, 'undefined function')
            canUseDeval = false;
        else
            rethrow(err);
        end
    end
end

xqRow = xq(:).';
yInterp = interp1(sol.x, sol.y.', xqRow, 'pchip', 'extrap').';

if isnumeric(xq) && isscalar(xq)
    yq = yInterp(:);
else
    yq = yInterp;
end
end

function emitRunLog(fmt, varargin)
txt = sprintf(fmt, varargin{:});
fprintf('%s', txt);
appendRunLog(txt);
end

function idx = fetchNextLogIndex()
idx = fetchAppCounter('pehf_logCounter');
end

function idx = beginLineLog()
idx = fetchAppCounter('pehf_lineLogCounter');
emitRunLog('[L%dS] -----------------------\n', idx);
end

function endLineLog(idx)
emitRunLog('[L%dE] -----------------------\n', idx);
end

function logContinuousRangeSummary(ranges, sweepName, domainLabel, statusLabel)
if nargin < 2 || isempty(sweepName)
    sweepName = '';
end
if nargin < 3 || isempty(domainLabel)
    domainLabel = 'parameter';
end
sweepLabel = sanitizeSectionTitle(sweepName);
paramLabel = sanitizeAxisLabelLocal(domainLabel);
statusLabel = lower(strtrim(string(statusLabel)));
if strlength(statusLabel) == 0
    statusLabel = 'status';
end
if isempty(ranges)
    emitRunLog('Continuous %s ranges (%s): none\n', statusLabel, sweepLabel);
    return
end
emitRunLog('Continuous %s ranges (%s):\n', statusLabel, sweepLabel);

for idx = 1:numel(ranges)
    range = ranges(idx);
    valueTxt = buildRangeText(range, paramLabel);
    emitRunLog('  Branch %d | attempts %d-%d | %s\n', ...
        range.branchIdx, range.startAttempt, range.endAttempt, valueTxt);
end
end

function label = sanitizeSectionTitle(name)
txt = char(string(name));
if strlength(string(txt)) == 0
    label = 'sweep';
else
    label = txt;
end
end

function logPlannedPointCount(sweepName, values, nBranches)
if nargin < 3
    return
end
if isempty(values)
    sweepCount = 0;
elseif numel(values) == 1 && isnan(values(1))
    sweepCount = 1;
else
    sweepCount = numel(values);
end
totalPoints = sweepCount * max(1, nBranches);
sweepLabel = sanitizeSectionTitle(sweepName);
emitRunLog('Planned points (%s): sweeps=%d, branches=%d, total=%d\n', ...
    sweepLabel, sweepCount, nBranches, totalPoints);
end

function txt = buildRangeText(range, paramLabel)
startVal = range.startValue;
endVal = range.endValue;
haveNumeric = isfinite(startVal) && isfinite(endVal);
if haveNumeric
    if abs(startVal - endVal) < eps(max(abs(startVal), abs(endVal))) * 10
        txt = sprintf('%s=%s', paramLabel, formatNumericToken(startVal));
    else
        txt = sprintf('%s=%s to %s', paramLabel, formatNumericToken(startVal), formatNumericToken(endVal));
    end
else
    startLabel = strtrim(char(range.startLabel));
    endLabel = strtrim(char(range.endLabel));
    if strlength(string(startLabel)) == 0 && strlength(string(endLabel)) == 0
        txt = 'values unavailable';
        return
    end
    if strcmp(startLabel, endLabel)
        txt = startLabel;
    else
        if strlength(string(startLabel)) == 0
            startLabel = 'start';
        end
        if strlength(string(endLabel)) == 0
            endLabel = 'end';
        end
        txt = sprintf('%s -> %s', startLabel, endLabel);
    end
end
end

function tokens = buildAttemptTokens(methodMeta, diagInfo)
tokens = {};
if nargin < 1 || isempty(methodMeta)
    methodMeta = struct();
end
if nargin < 2 || isempty(diagInfo)
    diagInfo = struct();
end
if isfield(methodMeta,'groupParamName') && ~isempty(methodMeta.groupParamName)
    val = NaN;
    if isfield(methodMeta,'groupParamValue')
        val = methodMeta.groupParamValue;
    end
    tokens{end+1} = formatParamToken(methodMeta.groupParamName, val); %#ok<AGROW>
end
paramLabel = 'x';
if isfield(methodMeta,'paramLabel') && ~isempty(methodMeta.paramLabel)
    paramLabel = methodMeta.paramLabel;
end
primaryVal = NaN;
if isfield(methodMeta,'domainDisplayValue') && isfinite(methodMeta.domainDisplayValue)
    primaryVal = methodMeta.domainDisplayValue;
elseif isfield(diagInfo,'domainValue') && isfinite(diagInfo.domainValue)
    primaryVal = diagInfo.domainValue;
end
if isfinite(primaryVal)
    tokens{end+1} = formatParamToken(paramLabel, primaryVal); %#ok<AGROW>
end
if isfield(methodMeta,'meshLabel') && ~isempty(methodMeta.meshLabel) ...
        && isfield(diagInfo,'domainRange') && all(isfinite(diagInfo.domainRange))
    tokens{end+1} = formatRangeToken(methodMeta.meshLabel, diagInfo.domainRange); %#ok<AGROW>
end
tokens = tokens(~cellfun('isempty', tokens));
end

function token = formatParamToken(label, value)
lbl = sanitizeAxisLabelLocal(label);
if isfinite(value)
    token = sprintf('%s=%s', lbl, formatNumericToken(value));
else
    token = '';
end
end

function token = formatRangeToken(label, rangeVals)
if numel(rangeVals) ~= 2 || any(~isfinite(rangeVals))
    token = '';
    return
end
lbl = sanitizeAxisLabelLocal(label);
token = sprintf('%s=[%s,%s]', lbl, formatNumericToken(rangeVals(1)), formatNumericToken(rangeVals(2)));
end

function token = formatNumericToken(value)
token = numericFormat('token', value);
end

function txt = formatNumericArray(values)
txt = numericFormat('array', values);
end

function outliers = detectSuccessOutliers(attemptLog, outlierSigma)
outliers = struct('branchIdx',{},'attemptIndex',{},'value',{},'label',{},'primaryMetric',{},'secondaryMetric',{}, ...
    'metricField',{},'median',{},'scale',{},'deviation',{},'normalizedDeviation',{});
if nargin == 0 || isempty(attemptLog)
    return
end
if nargin < 2 || isempty(outlierSigma) || ~isfinite(outlierSigma) || outlierSigma <= 0
    outlierSigma = 3;
end
branchCount = 0;
for idx = 1:numel(attemptLog)
    branchCount = max(branchCount, numel(attemptLog(idx).branchStatus));
end
for branchIdx = 1:branchCount
    values = [];
    secValues = [];
    attemptIdxList = [];
    sweepVals = [];
    labels = strings(0);
    for attemptIdx = 1:numel(attemptLog)
        entry = attemptLog(attemptIdx);
        statusVal = '';
        if branchIdx <= numel(entry.branchStatus) && ~isempty(entry.branchStatus{branchIdx})
            statusVal = lower(strtrim(string(entry.branchStatus{branchIdx})));
        end
        if ~strcmp(statusVal, 'success')
            continue
        end
        primaryVal = extractBranchValue(entry.solutionValues, branchIdx);
        secondaryVal = extractBranchValue(entry.secondaryValues, branchIdx);
        if ~isfinite(primaryVal) && ~isfinite(secondaryVal)
            continue
        end
        values(end+1) = primaryVal; %#ok<AGROW>
        secValues(end+1) = secondaryVal; %#ok<AGROW>
        attemptIdxList(end+1) = attemptIdx; %#ok<AGROW>
        sweepVals(end+1) = normalizeSweepValue(entry.value); %#ok<AGROW>
        labels(end+1) = string(entry.label); %#ok<AGROW>
    end
    [outliers, ~] = appendMetricOutliers(outliers, values, secValues, attemptIdxList, sweepVals, labels, branchIdx, outlierSigma, true);
    [outliers, ~] = appendMetricOutliers(outliers, values, secValues, attemptIdxList, sweepVals, labels, branchIdx, outlierSigma, false);
end
end

function [outliers, anyAdded] = appendMetricOutliers(outliers, values, secValues, attemptIdxList, sweepVals, labels, branchIdx, outlierSigma, usePrimary)
anyAdded = false;
if usePrimary
    metricVals = values;
else
    metricVals = secValues;
end

function field = ternaryMetricField(usePrimary)
if usePrimary
    field = 'Cf_star';
else
    field = 'Nu_star';
end
end
finiteMask = isfinite(metricVals);
metricVals = metricVals(finiteMask);
if numel(metricVals) < 3
    return
end
medVal = median(metricVals);
absDev = abs(metricVals - medVal);
madVal = median(absDev);
scale = 1.4826 * madVal;
if scale < 1e-9
    scale = std(metricVals, 1);
end
if scale < 1e-9
    return
end
threshold = outlierSigma * scale;
idxList = find(finiteMask);
for k = 1:numel(idxList)
    idx = idxList(k);
    deviation = abs(metricVals(k) - medVal);
    if deviation > threshold
        entryStruct = struct( ...
            'branchIdx', branchIdx, ...
            'attemptIndex', attemptIdxList(idx), ...
            'value', sweepVals(idx), ...
            'label', labels(idx), ...
            'primaryMetric', values(idx), ...
            'secondaryMetric', secValues(idx), ...
            'metricField', ternaryMetricField(usePrimary), ...
            'median', medVal, ...
            'scale', scale, ...
            'deviation', deviation, ...
            'normalizedDeviation', deviation / scale);
        outliers(end+1) = entryStruct; %#ok<AGROW>
        anyAdded = true;
    end
end
end

function val = extractBranchValue(arrayValues, branchIdx)
val = NaN;
if isempty(arrayValues)
    return
end
if branchIdx > numel(arrayValues)
    return
end
candidate = arrayValues(branchIdx);
if isnumeric(candidate) && isscalar(candidate)
    val = candidate;
end
end

function logOutlierSummary(outliers, sweepName, domainLabel)
sweepLabel = sanitizeSectionTitle(sweepName);
paramLabel = sanitizeAxisLabelLocal(domainLabel);
if isempty(outliers)
    emitRunLog('Success outliers (%s): none\n', sweepLabel);
    emitRunLog('=======================\n');
    return
end
emitRunLog('Success outliers (%s):\n', sweepLabel);
for idx = 1:numel(outliers)
    entry = outliers(idx);
    labelTxt = strtrim(char(entry.label));
    if strlength(string(labelTxt)) == 0
        labelTxt = sprintf('%s=%s', paramLabel, formatNumericToken(entry.value));
    end
    metricLabel = 'Cf';
    metricVal = entry.primaryMetric;
    if isfield(entry,'metricField') && strcmpi(entry.metricField, 'Nu_star')
        metricLabel = 'Nu';
        metricVal = entry.secondaryMetric;
    end
    emitRunLog('  Branch %d | attempt %d | %s | %s=%s | deviation=%s (%sx)\n', ...
        entry.branchIdx, entry.attemptIndex, labelTxt, metricLabel, formatNumericToken(metricVal), ...
        formatNumericToken(entry.deviation), formatNumericToken(entry.normalizedDeviation));
end

emitRunLog('=======================\n');
end

function val = normalizeSweepValue(value)
if isnumeric(value) && isscalar(value)
    val = value;
else
    val = NaN;
end
end

function fallback = buildMetricFallbacks(attemptLog)
fallback = struct('entries', {});
if isempty(attemptLog)
    return
end
maxBranches = 0;
for idx = 1:numel(attemptLog)
    maxBranches = max(maxBranches, numel(attemptLog(idx).branchStatus));
end
fallback = repmat(struct('entries', struct('value',{},'label',{},'metrics',struct())), maxBranches, 1);
for idx = 1:numel(attemptLog)
    entry = attemptLog(idx);
    sweepVal = entry.value;
    sweepLabel = entry.label;
    for b = 1:maxBranches
        statusVal = '';
        if b <= numel(entry.branchStatus) && ~isempty(entry.branchStatus{b})
            statusVal = lower(string(entry.branchStatus{b}));
        end
        if strcmp(statusVal, 'success')
            metrics = struct('Cf_star', NaN, 'Nu_star', NaN);
            if b <= numel(entry.solutionValues)
                metrics.Cf_star = entry.solutionValues(b);
            end
            if b <= numel(entry.secondaryValues)
                metrics.Nu_star = entry.secondaryValues(b);
            end
            fallback(b).entries(end+1) = struct( ... %#ok<AGROW>
                'value', sweepVal, ...
                'label', sweepLabel, ...
                'metrics', metrics);
        end
    end
end
end

function [metricVal, labelTxt] = fetchFallbackMetric(fallback, branchIdx, sweepValue, metricField)
metricVal = NaN;
labelTxt = '';
if nargin < 4
    metricField = 'Cf_star';
end
if branchIdx < 1 || branchIdx > numel(fallback)
    return
end
entries = fallback(branchIdx).entries;
if isempty(entries)
    return
end
tol = max(1e-8, 1e-8 * max(1, abs(sweepValue)));
for k = 1:numel(entries)
    entry = entries(k);
    if (isnan(sweepValue) && isnan(entry.value)) || abs(entry.value - sweepValue) <= tol
        if isfield(entry.metrics, metricField)
            metricVal = entry.metrics.(metricField);
        end
        labelTxt = entry.label;
        return
    end
end
end

function logMetricFallbackUsage(figCfg, branchIdx, sweepValue, sweepLabel)
paramLabel = sanitizeAxisLabelLocal(getfieldWithDefault(figCfg,'xlabel','parameter'));
labelTxt = sweepLabel;
if strlength(string(labelTxt)) == 0
    labelTxt = sprintf('%s=%s', paramLabel, formatNumericToken(sweepValue));
end
emitRunLog('Reused cached metrics for %s | branch %d due to missing results.\n', strtrim(char(labelTxt)), branchIdx);
end

function appendRunLog(txt)
try
    if ~isappdata(0,'pehf_logPath')
        return
    end
    logPath = getappdata(0,'pehf_logPath');
    if isempty(logPath)
        return
    end
    fid = fopen(logPath,'a');
    if fid < 0
        return
    end
    fprintf(fid,'%s', txt);
    fclose(fid);
catch
    % Logging failures should not affect solver runs.
end
end

function idx = fetchAppCounter(name)
idx = 1;
try
    if isappdata(0,name)
        idx = getappdata(0,name);
        setappdata(0,name, idx + 1);
    else
        setappdata(0,name, idx + 1);
    end
catch
    idx = 1;
end
end
