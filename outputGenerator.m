function outputInfo = outputGenerator(problemName, runOptions)
%OUTPUTGENERATOR Headless batch runner reusing the profile pipeline.

if nargin < 1 || strlength(string(problemName)) == 0
    problemName = 'natasha';
end
if nargin < 2 || isempty(runOptions)
    runOptions = struct();
end

cfg = factory('loadProblemValues', problemName);
model = factory('loadProblemModel', problemName);
method = factory('loadProblemMethod', problemName, cfg);

summary = factory('runProblem', cfg, model, method, struct('persistFiguresAsImages', true));

targetSpec = resolveTargetSweep(summary.sweeps, cfg, runOptions);
csvPath = writeSweepCsv(targetSpec, cfg, summary.outputDir, summary.figurePayloads);

outputInfo = struct( ...
    'profile', summary.problem, ...
    'outputDir', summary.outputDir, ...
    'sweepCsv', csvPath, ...
    'figureArtifacts', summary.figurePayloads);

fprintf('Problem "%s" artifacts stored in %s\n', summary.problem, summary.outputDir);
fprintf('Sweep CSV (%s): %s\n', targetSpec.displayName, csvPath);
end

function targetSpec = resolveTargetSweep(sweeps, cfg, opts)
if isempty(sweeps)
    error('outputGenerator:noSweep','No sweeps were executed for this profile.');
end
targetName = '';
if isfield(opts,'sweepName') && strlength(string(opts.sweepName)) > 0
    targetName = lower(string(opts.sweepName));
end
idx = 1;
if strlength(targetName) > 0
    names = string({sweeps.name});
    matches = find(lower(strtrim(names)) == targetName, 1, 'first');
    if ~isempty(matches)
        idx = matches;
    end
end
sweepEntry = sweeps(idx);
branchIdx = 1;
if isfield(opts,'branchIdx') && ~isempty(opts.branchIdx)
    branchIdx = opts.branchIdx;
end
metricField = 'Cf_star';
if isfield(opts,'metricField') && strlength(string(opts.metricField)) > 0
    metricField = char(opts.metricField);
end
groupValue = [];
if isfield(opts,'groupValue')
    groupValue = opts.groupValue;
elseif isfield(sweepEntry,'groupParamName')
    paramScope = 'p';
    if isfield(sweepEntry,'groupParamScope') && strlength(string(sweepEntry.groupParamScope)) > 0
        paramScope = lower(char(string(sweepEntry.groupParamScope)));
    end
    if strcmp(paramScope, 'm') && isfield(cfg,'m') && isfield(cfg.m, sweepEntry.groupParamName)
        groupValue = cfg.m.(sweepEntry.groupParamName);
    elseif isfield(cfg,'p') && isfield(cfg.p, sweepEntry.groupParamName)
        groupValue = cfg.p.(sweepEntry.groupParamName);
    end
end
selection = selectGroupEntry(sweepEntry, groupValue);
targetSpec = struct();
targetSpec.results = selection.results;
targetSpec.meta = selection.meta;
targetSpec.fileSuffix = selection.suffix;
targetSpec.displayName = string(sweepEntry.name);
targetSpec.branchIdx = branchIdx;
targetSpec.metricField = metricField;
targetSpec.groupLabel = selection.groupLabel;
targetSpec.sweepNameSlug = sanitizeFileToken(targetSpec.displayName);
end

function csvPath = writeSweepCsv(targetSpec, cfg, outDir, figurePayloads)
domainMinList = cfg.domainMinList;
domainMaxList = cfg.domainMaxList;
stepSizeList = cfg.domainGridSizeList;
results = targetSpec.results;
sweepMeta = targetSpec.meta;
suffix = targetSpec.fileSuffix;
plotStatusLookup = buildPlotStatusLookupOG(figurePayloads);
primaryCfg = findMetricConfigOG(figurePayloads, 'Cf_star');
secondaryCfg = findMetricConfigOG(figurePayloads, 'Nu_star');
successOutliers = [];
if isfield(sweepMeta, 'successOutliers') && ~isempty(sweepMeta.successOutliers)
    successOutliers = sweepMeta.successOutliers;
end

xVals = sweepMeta.sweepValues(:);
if isempty(xVals) || (numel(xVals) == 1 && isnan(xVals))
    error('outputGenerator:noSweepValues','No sweep values returned from solver.');
end

branchIdx = targetSpec.branchIdx;
yVals = extractBranchMetric(results, branchIdx, targetSpec.metricField, numel(xVals));
stepSpacing = computeStepSpacing(domainMinList, domainMaxList, stepSizeList, branchIdx);
stepSizes = repmat(stepSpacing, numel(xVals), 1);

successMask = isfinite(yVals);
xDisplay = arrayfun(@(xVal, isFail) formatXValue(xVal, isFail), xVals, ~successMask, 'UniformOutput', false);
statusValues = repmat("success", numel(xVals), 1);
statusValues(~successMask) = "failed";

csvPath = fullfile(outDir, sprintf('%s_branch%d_results%s.csv', targetSpec.sweepNameSlug, branchIdx, suffix));
fid = fopen(csvPath, 'w');
if fid == -1
    error('outputGenerator:writeSweepCsvFailed', 'Unable to open %s for writing.', csvPath);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
writeCsvRow(fid, {'domain_value','codomain_value','step_size','status'});
for idx = 1:numel(xDisplay)
    writeCsvRow(fid, {xDisplay{idx}, yVals(idx), stepSizes(idx), statusValues(idx)});
end
clear cleanupObj
if isfield(sweepMeta,'attemptLog')
    appendAttemptLogSection(csvPath, sweepMeta.attemptLog, branchIdx, plotStatusLookup, targetSpec.groupLabel, primaryCfg, secondaryCfg, successOutliers);
end
if isfield(sweepMeta,'continuousFailRanges')
    appendContinuousRangeSection(csvPath, sweepMeta.continuousFailRanges, branchIdx, 'fail');
end
if isfield(sweepMeta,'continuousSuccessRanges')
    appendContinuousRangeSection(csvPath, sweepMeta.continuousSuccessRanges, branchIdx, 'success');
end
if isfield(sweepMeta,'successOutliers')
    appendSuccessOutlierSection(csvPath, sweepMeta.successOutliers);
end
end

function yVals = extractBranchMetric(results, branchIdx, fieldName, numSweeps)
if isvector(results) && numel(results) ~= numSweeps
    numSweeps = numel(results);
end
yVals = nan(numSweeps, 1);
for k = 1:numSweeps
    entry = results(branchIdx, k);
    if isempty(entry.sol)
        continue
    end
    if isfield(entry, fieldName)
        yVals(k) = entry.(fieldName);
    end
end
end

function spacing = computeStepSpacing(domainMinList, domainMaxList, stepSizeList, branchIdx)
rangeStart = pickValue(domainMinList, branchIdx);
rangeEnd = pickValue(domainMaxList, branchIdx);
numPoints = max(2, round(pickValue(stepSizeList, branchIdx)));
spacing = (rangeEnd - rangeStart) / (numPoints - 1);
end

function val = pickValue(list, idx)
if isscalar(list)
    val = list;
else
    val = list(min(idx, numel(list)));
end
end

function outStr = formatXValue(xVal, isFail)
token = formatNumericToken(xVal);
if isFail
    outStr = sprintf('RED(%s)', token);
else
    outStr = token;
end
end

function appendAttemptLogSection(csvPath, attemptLog, branchIdx, plotStatusLookup, groupLabel, primaryCfg, secondaryCfg, successOutliers)
if nargin < 3
    branchIdx = 1;
end
if nargin < 4
    plotStatusLookup = struct();
end
if nargin < 5
    groupLabel = '';
end
if nargin < 6
    primaryCfg = struct();
end
if nargin < 7
    secondaryCfg = struct();
end
if nargin < 8
    successOutliers = [];
end

function appendContinuousRangeSection(csvPath, ranges, branchIdx, statusKey)
if nargin < 3
    branchIdx = 1;
end

function appendSuccessOutlierSection(csvPath, outliers)
if nargin < 2 || isempty(outliers)
    return
end
fid = fopen(csvPath, 'a');
if fid == -1
    warning('outputGenerator:appendOutlierFailed', 'Unable to append outlier section to %s.', csvPath);
    return
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
writeCsvRow(fid, {''});
writeCsvRow(fid, {'section','success_outliers'});
header = {'attempt_index','branch','sweep_value','sweep_label','metric_field','primary_metric','secondary_metric','deviation','normalized_deviation','median','scale'};
writeCsvRow(fid, header);
for idx = 1:numel(outliers)
    entry = outliers(idx);
    sweepLabel = formatRangeLabelOG(entry.label);
    metricField = getfieldWithDefault(entry, 'metricField', '');
    row = {entry.attemptIndex, entry.branchIdx, entry.value, sweepLabel, metricField, entry.primaryMetric, entry.secondaryMetric, entry.deviation, entry.normalizedDeviation, entry.median, entry.scale};
    writeCsvRow(fid, row);
end
writeCsvRow(fid, {''});
end
if nargin < 4 || strlength(string(statusKey)) == 0
    statusKey = 'status';
end
if isempty(ranges)
    return
end
fid = fopen(csvPath, 'a');
if fid == -1
    warning('outputGenerator:appendFailRangeFailed', 'Unable to append ranges to %s.', csvPath);
    return
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
writeCsvRow(fid, {''});
sectionKey = sprintf('continuous_%s_ranges_branch_%d', lower(strtrim(string(statusKey))), branchIdx);
writeCsvRow(fid, {'section', sectionKey});
header = {'branch_index','start_attempt','end_attempt','count','start_value','end_value','start_label','end_label'};
writeCsvRow(fid, header);
for idx = 1:numel(ranges)
    entry = ranges(idx);
    startLabel = formatRangeLabelOG(entry.startLabel);
    endLabel = formatRangeLabelOG(entry.endLabel);
    row = {entry.branchIdx, entry.startAttempt, entry.endAttempt, entry.count, ...
        entry.startValue, entry.endValue, startLabel, endLabel};
    writeCsvRow(fid, row);
end
writeCsvRow(fid, {''});
end
if isempty(attemptLog)
    return
end
fid = fopen(csvPath, 'a');
if fid == -1
    warning('outputGenerator:appendLogFailed', 'Unable to append attempt log to %s.', csvPath);
    return
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
writeCsvRow(fid, {''});
    writeCsvRow(fid, {'section', sprintf('attempt_log_branch_%d', branchIdx)});
    primaryLabel = buildMetricColumnLabelOG(attemptLog, branchIdx, 'primaryMetricLabels', 'primary_metric');
    secondaryLabel = buildMetricColumnLabelOG(attemptLog, branchIdx, 'secondaryMetricLabels', 'secondary_metric');
    header = {'attempt_index','sweep_value','sweep_label','all_branches_success','plot_b1_primary','plot_b2_primary','plot_b1_secondary','plot_b2_secondary', ...
        'plot_b1_primary_message','plot_b2_primary_message','plot_b1_secondary_message','plot_b2_secondary_message', ...
        'branches_deviation','branches_deviation_cf','branches_deviation_nu','status','configured_step','sweep_step','avg_mesh_step','max_residual','absolute_error','relative_error','error_message', ...
        'initial_guess_error','iterations', primaryLabel, secondaryLabel, 'primary_rate', 'secondary_rate'};
writeCsvRow(fid, header);
for idx = 1:numel(attemptLog)
    entry = attemptLog(idx);
    statusVal = safeCellFetch(entry.branchStatus, branchIdx, '');
    cfgStep = safeArrayValue(entry.stepSizes, branchIdx, NaN);
    avgStep = safeArrayValue(entry.avgMeshSteps, branchIdx, NaN);
    maxResidual = safeArrayValue(entry.maxResiduals, branchIdx, NaN);
    absError = safeArrayValue(getfieldWithDefault(entry, 'absoluteErrors', []), branchIdx, NaN);
    relError = safeArrayValue(getfieldWithDefault(entry, 'relativeErrors', []), branchIdx, NaN);
    errMsg = safeCellFetch(entry.errorMessages, branchIdx, '');
    if isSkipMessage(errMsg)
        statusVal = NaN;
    end
    solVal = safeArrayValue(entry.solutionValues, branchIdx, NaN);
    secVal = safeArrayValue(entry.secondaryValues, branchIdx, NaN);
    primaryRate = safeArrayValue(getfieldWithDefault(entry, 'primaryRates', []), branchIdx, NaN);
    secondaryRate = safeArrayValue(getfieldWithDefault(entry, 'secondaryRates', []), branchIdx, NaN);
    sweepStep = safeArrayValue(entry.sweepSteps, branchIdx, NaN);
    guessErr = safeArrayValue(entry.initialGuessErrors, branchIdx, NaN);
    iterCount = safeArrayValue(entry.iterations, branchIdx, NaN);
    if isfield(entry, 'allBranchesSuccess') && ~isempty(entry.allBranchesSuccess)
        allBranchesSuccess = logical(entry.allBranchesSuccess);
    else
        allBranchesSuccess = false;
    end
    [plotB1Primary, plotB2Primary, plotB1Secondary, plotB2Secondary] = resolvePlotStatusFlagsOG(plotStatusLookup, groupLabel, entry.value);
    msgB1Primary = buildPlotStatusReasonOG(plotStatusLookup, groupLabel, entry, idx, 1, 'Cf_star', plotB1Primary, primaryCfg, successOutliers);
    msgB2Primary = buildPlotStatusReasonOG(plotStatusLookup, groupLabel, entry, idx, 2, 'Cf_star', plotB2Primary, primaryCfg, successOutliers);
    msgB1Secondary = buildPlotStatusReasonOG(plotStatusLookup, groupLabel, entry, idx, 1, 'Nu_star', plotB1Secondary, secondaryCfg, successOutliers);
    msgB2Secondary = buildPlotStatusReasonOG(plotStatusLookup, groupLabel, entry, idx, 2, 'Nu_star', plotB2Secondary, secondaryCfg, successOutliers);
    row = {idx, entry.value, entry.label, allBranchesSuccess, plotB1Primary, plotB2Primary, plotB1Secondary, plotB2Secondary, ...
        msgB1Primary, msgB2Primary, msgB1Secondary, msgB2Secondary, ...
        getfieldWithDefault(entry, 'branchesDeviation', NaN), ...
        getfieldWithDefault(entry, 'branchesDeviationPrimary', NaN), getfieldWithDefault(entry, 'branchesDeviationSecondary', NaN), ...
        statusVal, cfgStep, sweepStep, avgStep, maxResidual, absError, relError, errMsg, guessErr, iterCount, solVal, secVal, primaryRate, secondaryRate};
    writeCsvRow(fid, row);
end
end

function value = safeArrayValue(arr, idx, defaultVal)
if isempty(arr) || idx < 1 || idx > numel(arr) || isempty(arr(idx))
    value = defaultVal;
else
    value = arr(idx);
end
end

function value = safeCellFetch(arr, idx, defaultVal)
if isempty(arr) || idx < 1 || idx > numel(arr) || isempty(arr{idx})
    value = defaultVal;
else
    value = arr{idx};
end
end

function tf = isSkipMessage(msg)
tf = false;
if strlength(string(msg)) == 0
    return
end
msgLower = lower(string(msg));
tf = contains(msgLower, "skipped remaining");
end

function writeCsvRow(fid, cells)
if isempty(cells)
    fprintf(fid, '\n');
    return
end

function txt = formatRangeLabelOG(label)
if isstring(label)
    txt = char(label);
elseif ischar(label)
    txt = label;
else
    txt = char(string(label));
end
end
tokens = cellfun(@formatCsvToken, cells, 'UniformOutput', false);
fprintf(fid, '%s\n', strjoin(tokens, ','));
end

function token = formatCsvToken(value)
if isstring(value)
    token = quoteString(char(value));
elseif ischar(value)
    token = quoteString(value);
elseif isempty(value)
    token = '';
elseif isnumeric(value)
    if numel(value) ~= 1
        token = '';
    elseif isnan(value)
        token = 'NaN';
    elseif isinf(value)
        if value > 0
            token = 'Inf';
        else
            token = '-Inf';
        end
    else
        token = formatNumericToken(value);
    end
else
    token = quoteString(char(string(value)));
end
end

function token = formatNumericToken(value)
token = numericFormat('token', value);
end

function out = quoteString(str)
if isempty(str)
    out = '""';
else
    out = ['"', strrep(str, '"', '""'), '"'];
end
end

function selection = selectGroupEntry(sweepEntry, desiredValue)
if ~isfield(sweepEntry,'groupSummaries') || isempty(sweepEntry.groupSummaries)
    selection = struct('results', sweepEntry.results, 'meta', sweepEntry.meta, 'suffix', '', 'groupLabel', '');
    return
end
groups = sweepEntry.groupSummaries;
values = [groups.value];
if isempty(desiredValue)
    if isfield(sweepEntry,'primaryGroupIndex') && ~isempty(sweepEntry.primaryGroupIndex)
        idx = sweepEntry.primaryGroupIndex;
    else
        idx = 1;
    end
else
    [~,idx] = min(abs(values - desiredValue));
end
idx = max(1, min(idx, numel(groups)));
target = groups(idx);
suffix = '';
if isfield(sweepEntry,'groupParamName') && strlength(string(sweepEntry.groupParamName)) > 0 && isfinite(target.value)
    suffix = sprintf('_%s%s', sanitizeFileToken(sweepEntry.groupParamName), formatGroupValue(target.value));
end
selection = struct('results', target.results, 'meta', target.meta, 'suffix', suffix, 'groupLabel', target.label);
end

function label = buildMetricColumnLabelOG(attemptLog, branchIdx, fieldName, defaultLabel)
if nargin < 4 || isempty(defaultLabel)
    defaultLabel = 'metric';
end

function token = sanitizeFileToken(value)
strVal = strtrim(char(string(value)));
if isempty(strVal)
    token = 'sweep';
    return
end
strVal = lower(strVal);
strVal = regexprep(strVal, '[^a-z0-9]+', '_');
strVal = regexprep(strVal, '_+', '_');
strVal = regexprep(strVal, '^_|_$', '');
if isempty(strVal)
    token = 'sweep';
else
    token = strVal;
end
end

function txt = formatGroupValue(val)
if isnumeric(val)
    txt = strrep(sprintf('%.3g', val), '.', 'p');
else
    txt = sanitizeFileToken(val);
end
end
raw = extractMetricLabelOG(attemptLog, branchIdx, fieldName);
if strlength(string(raw)) == 0
    raw = defaultLabel;
end
label = sanitizeMetricLabelOG(raw);
end

function label = extractMetricLabelOG(attemptLog, branchIdx, fieldName)
label = '';
if isempty(attemptLog)
    return
end
for idx = 1:numel(attemptLog)
    entry = attemptLog(idx);
    if ~isfield(entry, fieldName)
        continue
    end
    labels = entry.(fieldName);
    if isempty(labels) || branchIdx > numel(labels)
        continue
    end
    candidate = labels{branchIdx};
    if strlength(string(candidate)) > 0
        label = candidate;
        return
    end
end
end

function txt = sanitizeMetricLabelOG(labelIn)
text = lower(strtrim(string(labelIn)));
if strlength(text) == 0
    txt = 'metric';
    return
end
text = regexprep(text, '[^a-z0-9]+', '_');
text = regexprep(text, '_+', '_');
text = regexprep(text, '^_|_$', '');
if strlength(text) == 0
    txt = 'metric';
else
    txt = char(text);
end
end

function cfg = findMetricConfigOG(payloads, metricField)
cfg = struct();
if nargin < 1 || isempty(payloads) || strlength(string(metricField)) == 0
    return
end
payloads = payloads(:);
for idx = 1:numel(payloads)
    entry = payloads{idx};
    if isempty(entry) || ~isfield(entry,'mode') || ~strcmpi(string(entry.mode), 'metric')
        continue
    end
    if ~isfield(entry,'config') || ~isfield(entry.config,'metricField')
        continue
    end
    if strcmpi(string(entry.config.metricField), string(metricField))
        cfg = entry.config;
        return
    end
end
end

function reason = buildPlotStatusReasonOG(plotStatusLookup, groupLabel, attempt, attemptIdx, branchIdx, metricField, isPlotted, figCfg, successOutliers)
reason = '';
if isPlotted
    return
end
if nargin < 4
    attemptIdx = NaN;
end
if nargin < 5
    branchIdx = 1;
end
if nargin < 6 || strlength(string(metricField)) == 0
    metricField = 'Cf_star';
end
if nargin < 7
    isPlotted = false;
end
if nargin < 8
    figCfg = struct();
end
if nargin < 9
    successOutliers = [];
end
if isempty(plotStatusLookup)
    reason = 'no_plot_data';
    return
end
if isempty(findPlotStatusGroupOG(plotStatusLookup, normalizePlotGroupLabelOG(groupLabel)))
    reason = 'group_not_found';
    return
end
statusVal = '';
if isfield(attempt,'branchStatus') && branchIdx <= numel(attempt.branchStatus)
    statusVal = lower(strtrim(string(attempt.branchStatus{branchIdx})));
end
errMsg = '';
if isfield(attempt,'errorMessages') && branchIdx <= numel(attempt.errorMessages)
    errMsg = string(attempt.errorMessages{branchIdx});
end
if isSkipMessage(errMsg)
    reason = 'skipped';
    return
end
if strlength(statusVal) == 0 || ~strcmp(statusVal, 'success')
    reason = 'solver_failed';
    return
end
metricVal = pickAttemptMetricValueOG(attempt, branchIdx, metricField);
if ~isfinite(metricVal)
    reason = 'metric_missing';
    return
end
if shouldHideZeroDeviationFromCfgOG(figCfg)
    if isZeroDeviationAttemptOG(attempt, metricField, figCfg)
        reason = 'zero_deviation_filtered';
        return
    end
end
if isOutlierFilteredAttemptOG(attempt, attemptIdx, branchIdx, metricField, metricVal, successOutliers, figCfg)
    reason = 'outlier_filtered';
    return
end
reason = 'not_plotted';
end

function val = pickAttemptMetricValueOG(attempt, branchIdx, metricField)
val = NaN;
if strcmpi(metricField, 'Nu_star')
    val = safeArrayValue(getfieldWithDefault(attempt, 'secondaryValues', []), branchIdx, NaN);
else
    val = safeArrayValue(getfieldWithDefault(attempt, 'solutionValues', []), branchIdx, NaN);
end
end

function tf = isZeroDeviationAttemptOG(attempt, metricField, figCfg)
tf = false;
if nargin < 3
    figCfg = struct();
end
if strcmpi(metricField, 'Nu_star')
    vals = getfieldWithDefault(attempt, 'secondaryValues', []);
else
    vals = getfieldWithDefault(attempt, 'solutionValues', []);
end
vals = vals(isfinite(vals));
if numel(vals) < 2
    return
end
deviation = computeBranchesDeviationOG(vals);
tol = getZeroDeviationToleranceFromCfgOG(figCfg, vals);
if isfinite(deviation) && deviation <= tol
    tf = true;
end
end

function deviation = computeBranchesDeviationOG(values)
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

function tf = shouldHideZeroDeviationFromCfgOG(figCfg)
tf = false;
if nargin < 1 || ~isstruct(figCfg)
    return
end
if isfield(figCfg,'plotZeroDeviation') && ~isempty(figCfg.plotZeroDeviation)
    val = figCfg.plotZeroDeviation;
    if islogical(val) && isscalar(val)
        tf = ~val;
    elseif isnumeric(val) && isscalar(val)
        tf = (val == 0);
    end
end
end

function tol = getZeroDeviationToleranceFromCfgOG(figCfg, vals)
tol = 1e-6;
if nargin < 1 || ~isstruct(figCfg)
    return
end
if isfield(figCfg,'zeroDeviationTolerance') && ~isempty(figCfg.zeroDeviationTolerance)
    candidate = figCfg.zeroDeviationTolerance;
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && candidate >= 0
        tol = candidate;
    end
end
if nargin >= 2 && ~isempty(vals)
    tol = max(tol, tol * max(1, max(abs(vals))));
end
end

function tf = isOutlierFilteredAttemptOG(attempt, attemptIdx, branchIdx, metricField, metricVal, successOutliers, figCfg)
tf = false;
if isempty(successOutliers)
    return
end
if shouldPlotOutliersFromCfgOG(figCfg)
    return
end
filterMode = getOutlierFilterModeFromCfgOG(figCfg);
restrictToBranch = shouldRestrictOutliersForPlotOG(filterMode);
outlierXY = [];
for k = 1:numel(successOutliers)
    entry = successOutliers(k);
    if ~isfield(entry,'value') || ~isfield(entry,'branchIdx')
        continue
    end
    if restrictToBranch && isfinite(branchIdx) && entry.branchIdx ~= branchIdx
        continue
    end
    if strcmpi(metricField, 'Nu_star')
        if ~isfield(entry,'secondaryMetric')
            continue
        end
        yVal = entry.secondaryMetric;
    else
        if ~isfield(entry,'primaryMetric')
            continue
        end
        yVal = entry.primaryMetric;
    end
    if isfinite(entry.value) && isfinite(yVal)
        outlierXY(end+1,:) = [entry.value, yVal]; %#ok<AGROW>
    end
end
if isempty(outlierXY)
    return
end
tol = 1e-9 * max(1, max(abs([metricVal, attempt.value])));
dx = abs(outlierXY(:,1) - attempt.value);
dy = abs(outlierXY(:,2) - metricVal);
tf = any(dx <= tol & dy <= tol);
if tf
    return
end
if isfinite(attemptIdx)
    for k = 1:numel(successOutliers)
        entry = successOutliers(k);
        if ~isfield(entry,'attemptIndex') || ~isfield(entry,'branchIdx')
            continue
        end
        if entry.attemptIndex ~= attemptIdx
            continue
        end
        if restrictToBranch && entry.branchIdx ~= branchIdx
            continue
        end
        tf = true;
        return
    end
end
end

function tf = shouldPlotOutliersFromCfgOG(figCfg)
tf = false;
if nargin < 1 || ~isstruct(figCfg) || ~isfield(figCfg,'showOutliers')
    return
end
val = figCfg.showOutliers;
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

function mode = getOutlierFilterModeFromCfgOG(figCfg)
mode = 'global';
if nargin < 1 || ~isstruct(figCfg)
    return
end
if isfield(figCfg,'outlierFilter') && ~isempty(figCfg.outlierFilter)
    mode = char(string(figCfg.outlierFilter));
end
end

function tf = shouldRestrictOutliersForPlotOG(filterMode)
mode = lower(strtrim(char(string(filterMode))));
tf = any(strcmp(mode, {'branch','per-branch','per_branch','perbranch'}));
end

function lookup = buildPlotStatusLookupOG(payloads)
lookup = struct('label',{},'primary',{},'secondary',{});
if nargin < 1 || isempty(payloads)
    return
end
payloads = payloads(:);
for idx = 1:numel(payloads)
    payload = payloads{idx};
    if isempty(payload) || ~isfield(payload,'mode') || ~strcmpi(string(payload.mode), 'metric')
        continue
    end
    if ~isfield(payload,'config') || ~isfield(payload.config,'metricField')
        continue
    end
    metricField = char(string(payload.config.metricField));
    if ~(strcmpi(metricField, 'Cf_star') || strcmpi(metricField, 'Nu_star'))
        continue
    end
    lineSets = payload.lineSets;
    if isempty(lineSets)
        continue
    end
    for ln = 1:numel(lineSets)
        line = lineSets(ln);
        if ~isfield(line,'branchIdx') || ~isfield(line,'x') || ~isfield(line,'y')
            continue
        end
        branchIdx = line.branchIdx;
        groupLabel = '';
        if isfield(line,'groupLabel') && ~isempty(line.groupLabel)
            groupLabel = line.groupLabel;
        end
        groupLabel = normalizePlotGroupLabelOG(groupLabel);
        [lookup, groupIdx] = upsertPlotStatusGroupOG(lookup, groupLabel);
        if strcmpi(metricField, 'Cf_star')
            lookup(groupIdx).primary = appendPlotValuesOG(lookup(groupIdx).primary, branchIdx, line.x, line.y);
        else
            lookup(groupIdx).secondary = appendPlotValuesOG(lookup(groupIdx).secondary, branchIdx, line.x, line.y);
        end
    end
end
end

function [plotB1Primary, plotB2Primary, plotB1Secondary, plotB2Secondary] = resolvePlotStatusFlagsOG(lookup, label, sweepValue)
plotB1Primary = false;
plotB2Primary = false;
plotB1Secondary = false;
plotB2Secondary = false;
if isempty(lookup) || ~isfinite(sweepValue)
    return
end
groupLabel = normalizePlotGroupLabelOG(label);
groupIdx = findPlotStatusGroupOG(lookup, groupLabel);
if isempty(groupIdx)
    return
end
primary = lookup(groupIdx).primary;
secondary = lookup(groupIdx).secondary;
plotB1Primary = isValuePlottedOG(fetchPlotValuesOG(primary, 1), sweepValue);
plotB2Primary = isValuePlottedOG(fetchPlotValuesOG(primary, 2), sweepValue);
plotB1Secondary = isValuePlottedOG(fetchPlotValuesOG(secondary, 1), sweepValue);
plotB2Secondary = isValuePlottedOG(fetchPlotValuesOG(secondary, 2), sweepValue);
end

function labelOut = normalizePlotGroupLabelOG(labelIn)
if strlength(string(labelIn)) == 0
    labelOut = '';
else
    labelOut = formatLegendLabelOG(labelIn);
end
end

function idx = findPlotStatusGroupOG(lookup, groupLabel)
idx = [];
for k = 1:numel(lookup)
    if strcmp(lookup(k).label, groupLabel)
        idx = k;
        return
    end
end
end

function [lookup, idx] = upsertPlotStatusGroupOG(lookup, groupLabel)
idx = findPlotStatusGroupOG(lookup, groupLabel);
if ~isempty(idx)
    return
end
idx = numel(lookup) + 1;
lookup(idx).label = groupLabel;
lookup(idx).primary = {};
lookup(idx).secondary = {};
end

function values = fetchPlotValuesOG(cellVals, branchIdx)
values = [];
if isempty(cellVals) || branchIdx < 1 || branchIdx > numel(cellVals)
    return
end
values = cellVals{branchIdx};
end

function cellVals = appendPlotValuesOG(cellVals, branchIdx, xVals, yVals)
if nargin < 4 || isempty(xVals) || isempty(yVals) || branchIdx < 1 || ~isfinite(branchIdx)
    return
end
mask = isfinite(xVals) & isfinite(yVals);
if ~any(mask)
    return
end
values = xVals(mask);
values = values(:);
if isempty(cellVals)
    cellVals = {};
end
if branchIdx > numel(cellVals)
    cellVals{branchIdx} = [];
end
cellVals{branchIdx} = unique([cellVals{branchIdx}(:); values]);
end

function tf = isValuePlottedOG(values, target)
tf = false;
if isempty(values) || ~isfinite(target)
    return
end
values = values(:);
tol = 1e-8 * max(1, max(abs([values; target])));
tf = any(abs(values - target) <= tol);
end

function out = formatLegendLabelOG(label)
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
out = replaceNumericTokensOG(txt);
end

function out = replaceNumericTokensOG(txt)
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
