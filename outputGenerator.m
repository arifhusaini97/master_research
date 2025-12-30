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
csvPath = writeSweepCsv(targetSpec, cfg, summary.outputDir);

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
elseif isfield(sweepEntry,'groupParamName') && isfield(cfg,'p') && isfield(cfg.p, sweepEntry.groupParamName)
    groupValue = cfg.p.(sweepEntry.groupParamName);
end
selection = selectGroupEntry(sweepEntry, groupValue);
targetSpec = struct();
targetSpec.results = selection.results;
targetSpec.meta = selection.meta;
targetSpec.fileSuffix = selection.suffix;
targetSpec.displayName = string(sweepEntry.name);
targetSpec.branchIdx = branchIdx;
targetSpec.metricField = metricField;
targetSpec.sweepNameSlug = sanitizeFileToken(targetSpec.displayName);
end

function csvPath = writeSweepCsv(targetSpec, cfg, outDir)
domainMinList = cfg.domainMinList;
domainMaxList = cfg.domainMaxList;
stepSizeList = cfg.domainGridSizeList;
results = targetSpec.results;
sweepMeta = targetSpec.meta;
suffix = targetSpec.fileSuffix;

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

csvTbl = table( ...
    string(xDisplay), ...
    yVals, ...
    stepSizes, ...
    statusValues, ...
    'VariableNames', {'domain_value','codomain_value','step_size','status'});

csvPath = fullfile(outDir, sprintf('%s_branch%d_results%s.csv', targetSpec.sweepNameSlug, branchIdx, suffix));
writetable(csvTbl, csvPath);
if isfield(sweepMeta,'attemptLog')
    appendAttemptLogSection(csvPath, sweepMeta.attemptLog, branchIdx);
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
token = sprintf('%.5g', xVal);
if isFail
    outStr = sprintf('RED(%s)', token);
else
    outStr = token;
end
end

function appendAttemptLogSection(csvPath, attemptLog, branchIdx)
if nargin < 3
    branchIdx = 1;
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
header = {'attempt_index','branch','sweep_value','sweep_label','primary_metric','secondary_metric','deviation','normalized_deviation','median','scale'};
writeCsvRow(fid, header);
for idx = 1:numel(outliers)
    entry = outliers(idx);
    sweepLabel = formatRangeLabelOG(entry.label);
    row = {entry.attemptIndex, entry.branchIdx, entry.value, sweepLabel, entry.primaryMetric, entry.secondaryMetric, entry.deviation, entry.normalizedDeviation, entry.median, entry.scale};
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
    header = {'attempt_index','sweep_value','sweep_label','status','configured_step','sweep_step','avg_mesh_step','max_residual','error_message', ...
        'initial_guess_error','iterations', primaryLabel, secondaryLabel};
writeCsvRow(fid, header);
for idx = 1:numel(attemptLog)
    entry = attemptLog(idx);
    statusVal = safeCellFetch(entry.branchStatus, branchIdx, '');
    cfgStep = safeArrayValue(entry.stepSizes, branchIdx, NaN);
    avgStep = safeArrayValue(entry.avgMeshSteps, branchIdx, NaN);
    maxResidual = safeArrayValue(entry.maxResiduals, branchIdx, NaN);
    errMsg = safeCellFetch(entry.errorMessages, branchIdx, '');
    solVal = safeArrayValue(entry.solutionValues, branchIdx, NaN);
    secVal = safeArrayValue(entry.secondaryValues, branchIdx, NaN);
    sweepStep = safeArrayValue(entry.sweepSteps, branchIdx, NaN);
    guessErr = safeArrayValue(entry.initialGuessErrors, branchIdx, NaN);
    iterCount = safeArrayValue(entry.iterations, branchIdx, NaN);
    row = {idx, entry.value, entry.label, statusVal, cfgStep, sweepStep, avgStep, maxResidual, errMsg, guessErr, iterCount, solVal, secVal};
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
        token = sprintf('%.12g', value);
    end
else
    token = quoteString(char(string(value)));
end
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
    selection = struct('results', sweepEntry.results, 'meta', sweepEntry.meta, 'suffix', '');
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
selection = struct('results', target.results, 'meta', target.meta, 'suffix', suffix);
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
