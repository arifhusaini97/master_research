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

nBranches  = numel(guesses);
blankTemplate = blankResultStruct();
numSweeps  = numel(odeHandles);
results    = repmat(blankTemplate, nBranches, numSweeps);

prevSolutions = cell(1, nBranches);
solutionCache = cell(nBranches, numSweeps);
processed = false(1, numSweeps);
hasSuccessfulPoint = false;
attemptLog = emptyAttemptLog();
lineBlockIdx = beginLineLog();

while true
    if isempty(values) || isempty(processed)
        break;
    end
    remainingOrder = remainingSolveOrder(values, processed, sweepCfgOut);
    if isempty(remainingOrder)
        break;
    end
    sweepIdx = remainingOrder(1);

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
    if successThisPoint
        hasSuccessfulPoint = true;
    end
    processed(sweepIdx) = true;

    if ~successThisPoint
        if shouldRefine(sweepCfgOut) && hasSuccessfulPoint
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

failRanges = computeContinuousFailRanges(attemptLog, 'fail');
successRanges = computeContinuousFailRanges(attemptLog, 'success');
successOutliers = detectSuccessOutliers(attemptLog);
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
opts = struct('doPlot', true, 'sweep', [], 'figureConfigs', [], 'metricEvaluators', struct(), 'domainLabel', 'x');
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
    labelFcn = @(val) sprintf('%s=%.3g', name, val);
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
    nearestSol = selectNearestSolution(solutionCache(branchIdx,:), values, idx);
    if isempty(nearestSol) && ~isempty(prevSolutions{branchIdx})
        nearestSol = prevSolutions{branchIdx};
    end

    [sol, usedContinuation, diag] = attemptSolveWithGuess(odeCurrent, bcCurrent, nearestSol, baseGuessFcn, ...
        domainMin_i, domainMax_i, stepSize_i, branchMethodCfg, domainLabel);
    branchDiagnostics(branchIdx) = diag;
    if isempty(sol) && usedContinuation
        [sol, usedContinuation, diag] = attemptSolveWithGuess(odeCurrent, bcCurrent, [], baseGuessFcn, ...
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
    CfVal = evaluateMetricWithFallback(metricEvaluators, 'Cf_star', sol, CfVal);
    NuVal = evaluateMetricWithFallback(metricEvaluators, 'Nu_star', sol, NuVal);

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
        labels{k} = cfg.labelFcn(vals(k));
    end
else
    for k = 1:numel(vals)
        labels{k} = sprintf('%g', vals(k));
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
        lineSets = renderMetricFigure(cfg, results, sweepValues, colors, lineStyles, numSweeps, nBranches, metricFallbacks, doPlot);
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
    if ~isempty(legendEntries)
        legend(legendEntries,'Location','best'); box on;
    end
end
end

function lineSets = renderMetricFigure(cfg, results, sweepValues, colors, lineStyles, numSweeps, nBranches, metricFallbacks, doPlot)
lineSets = struct('branchIdx',{},'lineLabel',{},'status',{},'x',{},'y',{});
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
        'y', []);
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
    yVals = real(yVals);
    validMask = isfinite(yVals) & isfinite(sweepValues);
    if ~any(validMask)
        lineSets(end+1) = lineStruct; %#ok<AGROW>
        continue
    end
    lineStruct.status = 'success';
    lineStruct.x = sweepValues(validMask).';
    lineStruct.y = yVals(validMask).';
    lineSets(end+1) = lineStruct; %#ok<AGROW>

    if doPlot
        ls = lineStyles{mod(branch-1, numel(lineStyles))+1};
        if cfg.useSharedColor
            baseColor = colors(1,:);
        else
            baseColor = colors(mod(idx-1, size(colors,1))+1,:);
        end
        branchColor = branchColorVariant(baseColor, branch);
        plot(lineStruct.x, lineStruct.y, 'LineWidth', 1.5, 'Color', branchColor, 'LineStyle', ls);
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
    if ~isempty(legendEntries)
        legend(legendEntries,'Location','best'); box on;
    end
end
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

function entry = composeLegendEntry(sweepLabel, branchIdx)
if isempty(sweepLabel)
    entry = sprintf('branch %d', branchIdx);
else
    entry = sprintf('%s, branch %d', sweepLabel, branchIdx);
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
entry.errorMessages = repmat({''}, 1, nBranches);
entry.meshPoints = nan(1, nBranches);
entry.usedContinuation = false(1, nBranches);
entry.solutionValues = nan(1, nBranches);
entry.secondaryValues = nan(1, nBranches);
entry.stateVectors = repmat({[]}, 1, nBranches);
entry.sweepSteps = nan(1, nBranches);
entry.initialGuessErrors = nan(1, nBranches);
entry.iterations = nan(1, nBranches);
entry.consoleLogs = repmat({''}, 1, nBranches);
entry.primaryMetricLabels = repmat({''}, 1, nBranches);
entry.secondaryMetricLabels = repmat({''}, 1, nBranches);
entry.initialGuessMesh = repmat({[]}, 1, nBranches);
entry.initialGuessProfiles = repmat({[]}, 1, nBranches);

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
end

function logStruct = emptyAttemptLog()
logStruct = struct('value',{},'label',{},'branchStatus',{}, ...
    'stepSizes',{},'avgMeshSteps',{},'maxResiduals',{}, ...
    'errorMessages',{},'meshPoints',{},'usedContinuation',{}, ...
    'solutionValues',{},'secondaryValues',{}, ...
    'primaryMetricLabels',{},'secondaryMetricLabels',{}, ...
    'sweepSteps',{},'initialGuessErrors',{},'iterations',{}, ...
    'consoleLogs',{},'stateVectors',{},'initialGuessMesh',{},'initialGuessProfiles',{});
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

function value = evaluateMetricWithFallback(metricEvaluators, fieldName, sol, fallback)
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
    candidate = fn(sol);
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
    formattedVals = strjoin(arrayfun(@(v) sprintf('%.6g', v), vec, 'UniformOutput', false), ', ');
    emitRunLog('Result: [%s]\n', formattedVals);
end
if nargin < 2 || ~isfinite(cfValue)
    emitRunLog('localSkinFriction: N/A\n');
else
    emitRunLog('localSkinFriction: %.6g\n', cfValue);
end
if nargin < 3 || ~isfinite(nuValue)
    emitRunLog('nusseltNumber: N/A\n');
else
    emitRunLog('nusseltNumber: %.6g\n', nuValue);
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
                valueSuffix = sprintf(' (%s=%.5g)', paramLabelClean, methodMeta.domainDisplayValue);
            end
            paramWindowTxt = sprintf('%s window: [%g, %g]%s;', paramLabelClean, window(1), window(2), valueSuffix);
        elseif isfield(methodMeta,'domainDisplayValue') && isfinite(methodMeta.domainDisplayValue)
            paramWindowTxt = sprintf('%s=%.5g;', paramLabelClean, methodMeta.domainDisplayValue);
        end
    elseif isfield(methodMeta,'domainDisplayValue') && isfinite(methodMeta.domainDisplayValue)
        paramWindowTxt = sprintf('%s=%.5g;', paramLabelClean, methodMeta.domainDisplayValue);
    end
    if isfield(methodMeta,'meshLabel') && ~isempty(methodMeta.meshLabel)
        meshLabel = sanitizeAxisLabelLocal(methodMeta.meshLabel);
    end
    if isfield(methodMeta,'solutionRange') && numel(methodMeta.solutionRange) >= 2
        solveBounds = methodMeta.solutionRange;
        if all(isfinite(solveBounds))
            meshSolvedTxt = sprintf('%s mesh solved range [%g, %g];', meshLabel, solveBounds(1), solveBounds(end));
        end
    end
    if isfield(methodMeta,'groupParamName') && ~isempty(methodMeta.groupParamName)
        groupLabel = sanitizeAxisLabelLocal(methodMeta.groupParamName);
        if isfield(methodMeta,'groupParamValue') && isfinite(methodMeta.groupParamValue)
            groupTxt = sprintf('%s=%.5g;', groupLabel, methodMeta.groupParamValue);
        else
            groupTxt = sprintf('%s group;', groupLabel);
        end
    end
end

meshWindowTxt = sprintf('%s mesh window: [%g, %g];', meshLabel, domainMin, domainMax);

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
    txt = mat2str(value, 6);
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
    sprintf('The maximum residual is  %.3e.', residualValue);
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
if isstruct(stats) && isfield(stats,'maxres') && ~isempty(stats.maxres)
    maxRes = stats.maxres;
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

function txt = buildRangeText(range, paramLabel)
startVal = range.startValue;
endVal = range.endValue;
haveNumeric = isfinite(startVal) && isfinite(endVal);
if haveNumeric
    if abs(startVal - endVal) < eps(max(abs(startVal), abs(endVal))) * 10
        txt = sprintf('%s=%.5g', paramLabel, startVal);
    else
        txt = sprintf('%s=%.5g to %.5g', paramLabel, startVal, endVal);
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
    token = sprintf('%s=%.5g', lbl, value);
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
token = sprintf('%s=[%g,%g]', lbl, rangeVals(1), rangeVals(2));
end

function outliers = detectSuccessOutliers(attemptLog)
outliers = struct('branchIdx',{},'attemptIndex',{},'value',{},'label',{},'primaryMetric',{},'secondaryMetric',{},'median',{},'scale',{},'deviation',{},'normalizedDeviation',{});
if nargin == 0 || isempty(attemptLog)
    return
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
        if ~isfinite(primaryVal)
            continue
        end
        secondaryVal = extractBranchValue(entry.secondaryValues, branchIdx);
        values(end+1) = primaryVal; %#ok<AGROW>
        secValues(end+1) = secondaryVal; %#ok<AGROW>
        attemptIdxList(end+1) = attemptIdx; %#ok<AGROW>
        sweepVals(end+1) = normalizeSweepValue(entry.value); %#ok<AGROW>
        labels(end+1) = string(entry.label); %#ok<AGROW>
    end
    if numel(values) < 3
        continue
    end
    medVal = median(values);
    absDev = abs(values - medVal);
    madVal = median(absDev);
    scale = 1.4826 * madVal;
    if scale < 1e-9
        scale = std(values, 1);
    end
    if scale < 1e-9
        continue
    end
    threshold = 3 * scale;
    for idx = 1:numel(values)
        deviation = abs(values(idx) - medVal);
        if deviation > threshold
            entryStruct = struct( ...
                'branchIdx', branchIdx, ...
                'attemptIndex', attemptIdxList(idx), ...
                'value', sweepVals(idx), ...
                'label', labels(idx), ...
                'primaryMetric', values(idx), ...
                'secondaryMetric', secValues(idx), ...
                'median', medVal, ...
                'scale', scale, ...
                'deviation', deviation, ...
                'normalizedDeviation', deviation / scale);
            outliers(end+1) = entryStruct; %#ok<AGROW>
        end
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
        labelTxt = sprintf('%s=%.5g', paramLabel, entry.value);
    end
    emitRunLog('  Branch %d | attempt %d | %s | Cf=%.6g | deviation=%.3g (%.2fx)\n', ...
        entry.branchIdx, entry.attemptIndex, labelTxt, entry.primaryMetric, entry.deviation, entry.normalizedDeviation);
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
    labelTxt = sprintf('%s=%.5g', paramLabel, sweepValue);
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
