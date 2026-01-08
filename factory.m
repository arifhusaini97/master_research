function varargout = factory(action, varargin)
%FACTORY Central dispatcher for shared Powell-Eyring utilities.

if nargin == 0
    error('factory:missingAction','An action name is required.');
end
actionKey = lower(strtrim(string(action)));
switch actionKey
    case {'loadproblemvalues','values'}
        [varargout{1:nargout}] = loadProblemValues(varargin{:});
    case {'loadproblemmodel','model'}
        [varargout{1:nargout}] = loadProblemModel(varargin{:});
    case {'runproblem','figuremanager'}
        [varargout{1:nargout}] = runProblem(varargin{:});
    case {'loadproblemmethod','method'}
        [varargout{1:nargout}] = loadProblemMethod(varargin{:});
    case {'evalsolution','eval'}
        [varargout{1:nargout}] = evalSolution(varargin{:});
    case {'rangelabel','range'}
        [varargout{1:nargout}] = buildRangeLabel(varargin{:});
    otherwise
        error('factory:unknownAction','Unknown factory action "%s".', action);
end
end

function cfg = loadProblemValues(problemName)
if nargin < 1 || strlength(string(problemName)) == 0
    problemName = 'natasha';
end
problemKey = char(lower(strtrim(string(problemName))));
cfg = invokeProfileFunction(problemKey, 'values');
cfg = applyProblemConfig(cfg, problemKey);
cfg.problemName = problemKey;
cfg = ensureDomainDefaults(cfg);
end

function method = loadProblemMethod(problemName, cfg)
if nargin < 1 || strlength(string(problemName)) == 0
    problemName = 'natasha';
end
problemKey = char(lower(strtrim(string(problemName))));
if nargin < 2
    cfg = struct();
end
methodFn = invokeProfileFunction(problemKey, 'method', cfg);
if isa(methodFn, 'function_handle')
    method = methodFn(cfg);
else
    method = methodFn;
end
method.problemName = problemKey;
method = ensureMethodDefaults(method, cfg);
end

function model = loadProblemModel(problemName)
if nargin < 1 || strlength(string(problemName)) == 0
    problemName = 'natasha';
end
problemKey = char(lower(strtrim(string(problemName))));
model = invokeProfileFunction(problemKey, 'profile');
model.problemName = problemKey;
end

function summary = runProblem(cfg, model, method, varargin)
if nargin < 3 || isempty(method)
    method = struct();
end
method = ensureMethodDefaults(method, cfg);
opts = parseProblemRunOptions(varargin{:});
cfg = ensureDomainDefaults(cfg);
cfg = applyMethodDomainConfig(cfg, method);
applyNumericFormatOptions(cfg);
if ~isfield(cfg,'seedLibrary')
    cfg.seedLibrary = struct();
end
baseHandles = model.createBaseHandles(cfg);
metricEvaluators = model.metricEvaluators(cfg);
figureConfigs = model.figureConfigs(cfg, metricEvaluators);
sweeps = model.sweeps(cfg, figureConfigs, metricEvaluators, baseHandles);

timestampStr = datestr(now,'yymmdd_HHMM');
methodSlug = sanitizeRunToken(getfieldWithDefault(method,'solver','bvp'));
plannedSummaries = buildPlannedSweepSummaries(sweeps);
runLabel = deriveOutputLabel(cfg, plannedSummaries);
outputDir = fullfile(opts.outputRoot, sprintf('%s_%s_%s_%s', timestampStr, cfg.problemName, methodSlug, runLabel));
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
logPath = fullfile(outputDir, 'run.log');
setappdata(0,'pehf_logPath', logPath);
logCleanup = onCleanup(@clearRunLogPath);
fidLog = fopen(logPath,'w');
if fidLog >= 0
    fclose(fidLog);
end
setappdata(0,'pehf_logCounter', 1);
setappdata(0,'pehf_lineLogCounter', 1);
setappdata(0,'pehf_figureLogCounter', 1);

sweepSummaries = repmat(struct('name',"", 'results',[], 'meta',[]), 1, numel(sweeps));
figurePayloads = {};

for idx = 1:numel(sweeps)
    sweepDef = sweeps(idx);
    sweepOptions = sweepDef.options;
    if ~isfield(sweepOptions, 'figureConfigs') || isempty(sweepOptions.figureConfigs)
        sweepOptions.figureConfigs = figureConfigs;
    end
    if ~isfield(sweepOptions, 'metricEvaluators') || isempty(sweepOptions.metricEvaluators)
        sweepOptions.metricEvaluators = metricEvaluators;
    end
    if ~isfield(sweepOptions,'domainLabel') || strlength(string(sweepOptions.domainLabel)) == 0
        sweepOptions.domainLabel = cfg.domainLabel;
    end

    if isfield(sweepDef, 'groupValues') && ~isempty(sweepDef.groupValues)
        [groupSummaries, groupedPayloads, primaryIdx] = runGroupedSweep(cfg, sweepDef, sweepOptions, model, method);
        sweepSummaries(idx).name = sweepDef.name;
        sweepSummaries(idx).results = groupSummaries(primaryIdx).results;
        sweepSummaries(idx).meta = groupSummaries(primaryIdx).meta;
        sweepSummaries(idx).groupSummaries = groupSummaries;
        sweepSummaries(idx).primaryGroupIndex = primaryIdx;
        if ~isempty(groupedPayloads)
            figurePayloads = [figurePayloads; groupedPayloads(:)]; %#ok<AGROW>
        end
        continue;
    end

    [sweepResults, sweepMeta, payloads] = displayBaseFn( ...
        baseHandles.ode, baseHandles.bc, ...
        cfg.guesses, cfg.domainMinList, cfg.domainMaxList, cfg.domainGridSizeList, ...
        method, sweepOptions);

    sweepSummaries(idx).name = sweepDef.name;
    sweepSummaries(idx).results = sweepResults;
    sweepSummaries(idx).meta = sweepMeta;
    cfg = updateSeedLibrary(cfg, sweepSummaries(idx));

    if ~isempty(payloads)
        figurePayloads = [figurePayloads; payloads(:)]; %#ok<AGROW>
    end
end

finalLabel = deriveOutputLabel(cfg, sweepSummaries);
if strlength(string(finalLabel)) == 0
    finalLabel = runLabel;
end
if ~strcmp(finalLabel, runLabel)
    candidateDir = fullfile(opts.outputRoot, sprintf('%s_%s_%s_%s', timestampStr, cfg.problemName, methodSlug, finalLabel));
    if ~strcmp(candidateDir, outputDir)
        if exist(candidateDir,'dir')
            candidateDir = uniquifyOutputDir(candidateDir);
        end
        movefile(outputDir, candidateDir);
        outputDir = candidateDir;
        logPath = fullfile(outputDir, 'run.log');
        setappdata(0,'pehf_logPath', logPath);
    end
    runLabel = finalLabel;
end

inputSnapshot = persistInputArtifacts(cfg.problemName, outputDir);

if ~isempty(figurePayloads)
    figurePayloads = attachParamSnapshot(figurePayloads, cfg);
    persistFigureArtifacts(figurePayloads, outputDir, opts.persistFiguresAsImages);
end

summary = struct();
summary.problem = cfg.problemName;
summary.outputDir = outputDir;
summary.sweeps = sweepSummaries;
summary.figurePayloads = figurePayloads;
summary.inputFiles = inputSnapshot;
end

function entries = buildPlannedSweepSummaries(sweeps)
if nargin == 0 || isempty(sweeps)
    entries = struct('meta', struct());
    return
end
entries = repmat(struct('meta', struct()), numel(sweeps), 1);
for idx = 1:numel(sweeps)
    meta = struct();
    if isfield(sweeps(idx),'name')
        meta.sweepName = sweeps(idx).name;
    else
        meta.sweepName = sprintf('sweep%d', idx);
    end
    meta.sweepValues = [];
    if isfield(sweeps(idx),'options') && isfield(sweeps(idx).options,'sweep') ...
            && isfield(sweeps(idx).options.sweep,'values')
        meta.sweepValues = sweeps(idx).options.sweep.values;
    end
    entries(idx).meta = meta;
end
end

function clearRunLogPath()
if isappdata(0,'pehf_logPath')
    rmappdata(0,'pehf_logPath');
end
if isappdata(0,'pehf_logCounter')
    rmappdata(0,'pehf_logCounter');
end
if isappdata(0,'pehf_lineLogCounter')
    rmappdata(0,'pehf_lineLogCounter');
end
if isappdata(0,'pehf_figureLogCounter')
    rmappdata(0,'pehf_figureLogCounter');
end
end

function newDir = uniquifyOutputDir(baseDir)
newDir = baseDir;
counter = 1;
while exist(newDir,'dir')
    newDir = [baseDir, '_', num2str(counter)];
    counter = counter + 1;
end
end

function opts = parseProblemRunOptions(arg)
opts = struct('outputRoot', 'output', 'persistFiguresAsImages', false);
if nargin == 0 || isempty(arg)
    return
end
if isstruct(arg)
    if isfield(arg,'outputRoot') && ~isempty(arg.outputRoot)
        opts.outputRoot = arg.outputRoot;
    end
    if isfield(arg,'persistFiguresAsImages') && ~isempty(arg.persistFiguresAsImages)
        opts.persistFiguresAsImages = logical(arg.persistFiguresAsImages);
    end
end
end

function out = invokeProfileFunction(problemName, funcStem, varargin)
funcName = char(funcStem);
problemDir = fullfile(projectRoot(), 'problems', problemName);
if ~isfolder(problemDir)
    error('factory:missingProblemDir','Problem folder not found: %s', problemDir);
end
addpath(problemDir);
cleaner = onCleanup(@() rmpath(problemDir));
if exist(funcName, 'file') ~= 2
    error('factory:missingFunction','Problem "%s" lacks function "%s.m".', problemName, funcName);
end
out = feval(funcName, varargin{:});
end

function root = projectRoot()
root = fileparts(mfilename('fullpath'));
end

function info = persistInputArtifacts(problemName, outDir)
info = struct('files',[]);
if nargin < 2 || strlength(string(problemName)) == 0 || strlength(string(outDir)) == 0
    return
end
probDir = fullfile(projectRoot(), 'problems', char(problemName));
if ~isfolder(probDir)
    return
end
targetDir = fullfile(outDir, 'input');
if ~exist(targetDir, 'dir')
    mkdir(targetDir);
end
fileList = {'method','model','profile','values','config'};
records = struct('name',{},'path',{});
for idx = 1:numel(fileList)
    src = fullfile(probDir, [fileList{idx}, '.m']);
    if exist(src,'file') ~= 2
        continue
    end
    dst = fullfile(targetDir, [fileList{idx}, '.m']);
    copyfile(src, dst);
    records(end+1) = struct('name', fileList{idx}, 'path', dst); %#ok<AGROW>
end
info.files = records;
end

% Objective: Merge optional problem config into values.
% Purpose: Keep values.m pure while still providing runtime settings.
% SWOT: S-clean separation; W-extra lookup; O-pluggable configs; T-missing config files.
function cfg = applyProblemConfig(cfg, problemKey)
cfg = ensureStruct(cfg);
configPath = fullfile(projectRoot(), 'problems', problemKey, 'config.m');
if exist(configPath, 'file') ~= 2
    return
end
cfgConfig = invokeProfileFunction(problemKey, 'config');
if ~isstruct(cfgConfig)
    return
end
cfg = mergeStructs(cfg, cfgConfig);
end

function out = mergeStructs(base, override)
out = base;
fields = fieldnames(override);
for k = 1:numel(fields)
    name = fields{k};
    out.(name) = override.(name);
end
end

function cfg = ensureStruct(cfg)
if ~isstruct(cfg)
    cfg = struct();
end
end

function label = buildRangeLabel(values)
if nargin == 0 || isempty(values)
    label = 'range';
    return
end
finiteVals = values(isfinite(values));
if isempty(finiteVals)
    label = 'range';
    return
end
label = sprintf('%sv%s', formatBound(min(finiteVals)), formatBound(max(finiteVals)));
end

function boundStr = formatBound(val)
prefix = '';
if val < 0
    prefix = 'm';
    val = abs(val);
end
raw = sprintf('%.3f', val);
raw = regexprep(raw, '0+$', '');
raw = regexprep(raw, '\.$', '');
if isempty(raw)
    raw = '0';
end
boundStr = [prefix, strrep(raw, '.', 'p')];
end

function info = persistFigureArtifacts(payloads, outDir, closeAfter)
if nargin < 3
    closeAfter = false;
end
payloads = payloads(:);
payloads = payloads(~cellfun('isempty', payloads));
if isempty(payloads)
    info = struct('figureId',{},'figPath',{},'pngPath',{},'csvPath',{}); %#ok<NASGU>
    return
end
figIds = cellfun(@(p) p.figureId, payloads);
savedFigures = persistFigures(figIds, outDir, closeAfter);
info = repmat(struct('figureId',NaN,'figPath','','pngPath','', 'csvPath',''), numel(payloads), 1);
for idx = 1:numel(payloads)
    figId = payloads{idx}.figureId;
    match = savedFigures([savedFigures.figureId] == figId);
    csvPath = persistFigureCsv(payloads{idx}, outDir);
    record = struct('figureId', figId, 'figPath', '', 'pngPath', '', 'csvPath', csvPath);
    if ~isempty(match)
        record.figPath = match(1).figPath;
        record.pngPath = match(1).pngPath;
    end
    info(idx) = record;
end
end

function saved = persistFigures(figIds, outDir, closeAfter)
if nargin < 3
    closeAfter = false;
end
figIds = unique(figIds(:).');
figHandles = findall(0,'Type','figure');
saved = repmat(struct('figureId',NaN,'figPath','','pngPath',''), numel(figIds), 1);
for idx = 1:numel(figIds)
    figId = figIds(idx);
    handle = locateFigure(figHandles, figId);
    if isempty(handle) || ~ishandle(handle)
        continue
    end
    try
        typeName = get(handle, 'Type');
    catch
        continue
    end
    if ~strcmpi(typeName, 'figure')
        continue
    end
    figPath = fullfile(outDir, sprintf('figure%d.fig', figId));
    pngPath = fullfile(outDir, sprintf('figure%d.png', figId));
    try
        savefig(handle, figPath);
        exportgraphics(handle, pngPath, 'Resolution', 300);
    catch err %#ok<NASGU>
        % Skip handles that cannot be persisted and continue with the rest.
        if exist(figPath,'file')
            delete(figPath);
        end
        if exist(pngPath,'file')
            delete(pngPath);
        end
        continue
    end
    saved(idx) = struct('figureId', figId, 'figPath', figPath, 'pngPath', pngPath);
    if closeAfter
        close(handle);
    end
end
end

function handle = locateFigure(figHandles, targetId)
handle = [];
for idx = 1:numel(figHandles)
    candidate = figHandles(idx);
    if isempty(candidate) || ~ishandle(candidate)
        continue
    end
    try
        figNumber = get(candidate,'Number');
    catch
        continue
    end
    if isequal(figNumber, targetId)
        handle = candidate;
        return
    end
end
end

function csvPath = persistFigureCsv(figPayload, outDir)
figId = figPayload.figureId;
cfg = figPayload.config;
csvPath = fullfile(outDir, sprintf('figure%d.csv', figId));
fid = fopen(csvPath,'w');
if fid < 0
    error('factory:writeFigureCsv','Unable to open %s for writing.', csvPath);
end
cleaner = onCleanup(@() fclose(fid));

writeCsvRow(fid, {'section','meta'});
writeCsvRow(fid, {'key','value'});
writeCsvRow(fid, {'figure_id', figId});
writeCsvRow(fid, {'title', cfg.title});
writeCsvRow(fid, {'mode', lower(string(cfg.mode))});
if isfield(cfg,'xlabel') && ~isempty(cfg.xlabel)
    writeCsvRow(fid, {'xlabel', cfg.xlabel});
end
if isfield(cfg,'ylabel') && ~isempty(cfg.ylabel)
    writeCsvRow(fid, {'ylabel', cfg.ylabel});
end
window = resolveFigureWindow(cfg);
if ~isempty(window)
    writeCsvRow(fid, {'domain_window', sprintf('[%s,%s]', formatNumericToken(window(1)), formatNumericToken(window(end)))});
end
if isfield(cfg,'codomainWindow') && ~isempty(cfg.codomainWindow)
    writeCsvRow(fid, {'codomain_window', sprintf('[%s,%s]', formatNumericToken(cfg.codomainWindow(1)), formatNumericToken(cfg.codomainWindow(end)))});
end
writeCsvRow(fid, {''});
if isfield(figPayload,'paramSnapshot') && ~isempty(figPayload.paramSnapshot)
    appendParamSnapshot(fid, figPayload.paramSnapshot);
end

lineSets = figPayload.lineSets;
if isempty(lineSets)
    writeCsvRow(fid, {'section','line'});
    writeCsvRow(fid, {'line_label','(none)'});
    writeCsvRow(fid, {'status','no_data'});
    return
end

for idx = 1:numel(lineSets)
    line = lineSets(idx);
    writeCsvRow(fid, {'section','line'});
    writeCsvRow(fid, {'line_index', idx});
    writeCsvRow(fid, {'line_label', line.lineLabel});
    if isfield(line,'branchIdx') && ~isempty(line.branchIdx)
        writeCsvRow(fid, {'branch_idx', line.branchIdx});
    end
    if isfield(line,'sweepLabel') && ~isempty(line.sweepLabel)
        writeCsvRow(fid, {'sweep_label', line.sweepLabel});
    end
    if isfield(line,'sweepValue') && isfinite(line.sweepValue)
        writeCsvRow(fid, {'sweep_value', line.sweepValue});
    end
    if isfield(line,'groupLabel') && ~isempty(line.groupLabel)
        writeCsvRow(fid, {'group_label', line.groupLabel});
    end
    if isfield(line,'asymptoteX')
        if isempty(line.asymptoteX)
            writeCsvRow(fid, {'asymptote_x', ''});
        else
            writeCsvRow(fid, {'asymptote_x', formatNumericArray(line.asymptoteX)});
        end
    end
    writeCsvRow(fid, {'status', line.status});

    if isempty(line.x)
        writeCsvRow(fid, {''});
        continue
    end

    xVals = line.x(:);
    yVals = line.y(:);
    stepSize = computePointSteps(xVals);

    xHeader = sanitizeAxisLabel(getfieldWithDefault(cfg,'xlabel','x'));
    if strcmpi(figPayload.mode, 'profile')
        columnHeader = {'domain_value','codomain_value','step_size','status'};
    else
        columnHeader = {'domain_value','codomain_value','step_size','status'};
    end
    writeCsvRow(fid, columnHeader);
    for row = 1:numel(xVals)
        writeCsvRow(fid, {xVals(row), yVals(row), stepSize(row), line.status});
    end
    writeCsvRow(fid, {''});
end

if isfield(figPayload,'groupLogs') && ~isempty(figPayload.groupLogs)
    for idx = 1:numel(figPayload.groupLogs)
        entry = figPayload.groupLogs(idx);
        figureLogIdx = fetchNextFigureLogIndex();
        emitFactoryLog('[F%dS] ***********************\n', figureLogIdx);
        writeAttemptLogSection(fid, entry.attemptLog, figPayload.branchCount, entry.label, figPayload.config);
        if isfield(entry,'failRanges')
            writeContinuousRangeSection(fid, entry.failRanges, entry.label, 'fail');
        end
        if isfield(entry,'successRanges')
            writeContinuousRangeSection(fid, entry.successRanges, entry.label, 'success');
        end
        if isfield(entry,'successOutliers')
            writeSuccessOutlierSection(fid, entry.successOutliers, entry.label);
        end
        emitFactoryLog('[F%dE] ***********************\n', figureLogIdx);
    end
elseif isfield(figPayload,'attemptLog') && ~isempty(figPayload.attemptLog)
    figureLogIdx = fetchNextFigureLogIndex();
    emitFactoryLog('[F%dS] ***********************\n', figureLogIdx);
    writeAttemptLogSection(fid, figPayload.attemptLog, figPayload.branchCount, '', figPayload.config);
    if isfield(figPayload,'failRanges')
        writeContinuousRangeSection(fid, figPayload.failRanges, '', 'fail');
    end
    if isfield(figPayload,'successRanges')
        writeContinuousRangeSection(fid, figPayload.successRanges, '', 'success');
    end
    if isfield(figPayload,'successOutliers')
        writeSuccessOutlierSection(fid, figPayload.successOutliers, '');
    end
    emitFactoryLog('[F%dE] ***********************\n', figureLogIdx);
end

end

function payloads = attachParamSnapshot(payloads, cfg)
if isempty(payloads) || ~isstruct(cfg)
    return
end
snapshot = struct();
if isfield(cfg,'p')
    snapshot.p = cfg.p;
end
if isfield(cfg,'m')
    snapshot.m = cfg.m;
end
if isfield(cfg,'n')
    snapshot.n = cfg.n;
end
if isempty(fieldnames(snapshot))
    return
end
for idx = 1:numel(payloads)
    if isempty(payloads{idx})
        continue
    end
    payloads{idx}.paramSnapshot = snapshot;
end
end

function appendParamSnapshot(fid, snapshot)
writeCsvRow(fid, {'section','parameters'});
writeCsvRow(fid, {'index','name','value'});
idx = 1;
idx = appendParamGroup(fid, snapshot, 'p', idx);
idx = appendParamGroup(fid, snapshot, 'm', idx);
appendParamGroup(fid, snapshot, 'n', idx);
writeCsvRow(fid, {''});
end

function idx = appendParamGroup(fid, snapshot, groupName, idx)
if ~isfield(snapshot, groupName)
    return
end
group = snapshot.(groupName);
if ~isstruct(group)
    return
end
names = fieldnames(group);
for k = 1:numel(names)
    name = names{k};
    value = group.(name);
    label = sprintf('%s.%s', groupName, name);
    writeCsvRow(fid, {idx, label, value});
    idx = idx + 1;
end
end

function [groupSummaries, combinedPayloads, primaryIdx] = runGroupedSweep(cfg, sweepDef, sweepOptions, model, method)
paramName = sweepDef.groupParamName;
groupVals = sweepDef.groupValues(:).';
labelFcn = sweepDef.groupLabelFcn;
numGroups = numel(groupVals);
figConfigs = sweepOptions.figureConfigs;
groupSummaries = repmat(struct('label',"",'value',NaN,'results',[],'meta',[]), numGroups, 1);
groupPayloads = cell(numGroups, numel(figConfigs));

for g = 1:numGroups
    cfgGroup = cfg;
    cfgGroup.p.(paramName) = groupVals(g);
    cfgGroup = ensureDomainDefaults(cfgGroup);
    if ~isfield(cfgGroup,'domainLabel') || isempty(cfgGroup.domainLabel) || strcmp(cfgGroup.domainLabel,'x')
        cfgGroup.domainLabel = getfieldWithDefault(sweepOptions,'domainLabel', cfg.domainLabel);
    end
    cfgGroup = applySeedGuesses(cfg, cfgGroup, paramName, groupVals(g));
    handles = model.createBaseHandles(cfgGroup);
    metrics = model.metricEvaluators(cfgGroup);
    if isfield(sweepDef,'optionsBuilder') && ~isempty(sweepDef.optionsBuilder)
        options = sweepDef.optionsBuilder(cfgGroup, figConfigs, metrics);
    else
        options = sweepOptions;
        options.metricEvaluators = metrics;
        options.figureConfigs = figConfigs;
    end
    options.metricEvaluators = metrics;
    options.figureConfigs = figConfigs;
    if ~isfield(options,'domainLabel') || strlength(string(options.domainLabel)) == 0
        options.domainLabel = getfieldWithDefault(sweepOptions,'domainLabel', cfgGroup.domainLabel);
    end
    options.doPlot = false;

    groupMethod = method;
    groupMethod.groupParamName = paramName;
    groupMethod.groupParamValue = groupVals(g);

    [sweepResults, sweepMeta, payloads] = displayBaseFn( ...
        handles.ode, handles.bc, ...
        cfgGroup.guesses, cfgGroup.domainMinList, cfgGroup.domainMaxList, cfgGroup.domainGridSizeList, ...
        groupMethod, options);

    groupSummaries(g).label = labelFcn(groupVals(g));
    groupSummaries(g).value = groupVals(g);
    groupSummaries(g).results = sweepResults;
    groupSummaries(g).meta = sweepMeta;
    for fIdx = 1:numel(payloads)
        groupPayloads{g, fIdx} = payloads{fIdx};
    end
end

combinedPayloads = combineGroupFigurePayloads(figConfigs, groupPayloads, groupSummaries);
primaryVal = cfg.p.(paramName);
[~,primaryIdx] = min(abs(groupVals - primaryVal));
if isempty(primaryIdx) || isnan(primaryIdx)
    primaryIdx = 1;
end
end

function payloads = combineGroupFigurePayloads(figConfigs, groupPayloads, groupSummaries)
numFigs = numel(figConfigs);
numGroups = numel(groupSummaries);
payloads = cell(0,1);
groupLabels = arrayfun(@(g) char(g.label), groupSummaries, 'UniformOutput', false);
groupLabels = cellfun(@formatLegendLabel, groupLabels, 'UniformOutput', false);
for figIdx = 1:numFigs
    cfg = figConfigs(figIdx);
    combinedLines = struct('branchIdx',{},'lineLabel',{},'status',{},'x',{},'y',{},'groupLabel',{},'sweepValue',{},'sweepLabel',{},'asymptoteX',{});
    branchCount = 0;
    for g = 1:numGroups
        entry = groupPayloads{g, figIdx};
        if isempty(entry)
            continue
        end
        branchCount = max(branchCount, entry.branchCount);
        lines = entry.lineSets;
        for ln = 1:numel(lines)
            line = lines(ln);
            lineStruct = struct( ...
                'branchIdx', getLineField(line,'branchIdx',NaN), ...
                'lineLabel', sprintf('%s | %s', groupLabels{g}, getLineField(line,'lineLabel','')), ...
                'status', getLineField(line,'status','no_data'), ...
                'x', getLineField(line,'x',[]), ...
                'y', getLineField(line,'y',[]), ...
                'sweepLabel', getLineField(line,'sweepLabel',''), ...
                'sweepValue', getLineField(line,'sweepValue',NaN), ...
                'groupLabel', groupLabels{g}, ...
                'asymptoteX', getLineField(line,'asymptoteX',[]));
            combinedLines(end+1) = lineStruct; %#ok<AGROW>
        end
    end
    if isempty(combinedLines)
        continue
    end
    plotGroupedMetricFigure(cfg, combinedLines, groupLabels);
    payloads{end+1} = struct( ... %#ok<AGROW>
        'config', cfg, ...
        'mode', cfg.mode, ...
        'figureId', cfg.figureId, ...
        'lineSets', {combinedLines}, ...
        'groupLogs', {buildGroupLogs(groupSummaries)}, ...
        'branchCount', branchCount);
end
end

function logs = buildGroupLogs(groupSummaries)
logs = struct('label',{},'attemptLog',{},'failRanges',{},'successRanges',{},'successOutliers',{});
for idx = 1:numel(groupSummaries)
    logs(idx).label = groupSummaries(idx).label;
    logs(idx).attemptLog = groupSummaries(idx).meta.attemptLog;
    if isfield(groupSummaries(idx).meta,'continuousFailRanges')
        logs(idx).failRanges = groupSummaries(idx).meta.continuousFailRanges;
    else
        logs(idx).failRanges = struct('branchIdx',{},'startAttempt',{},'endAttempt',{},'count',{},'startValue',{},'endValue',{},'startLabel',{},'endLabel',{});
    end
    if isfield(groupSummaries(idx).meta,'continuousSuccessRanges')
        logs(idx).successRanges = groupSummaries(idx).meta.continuousSuccessRanges;
    else
        logs(idx).successRanges = struct('branchIdx',{},'startAttempt',{},'endAttempt',{},'count',{},'startValue',{},'endValue',{},'startLabel',{},'endLabel',{});
    end
    if isfield(groupSummaries(idx).meta,'successOutliers')
        logs(idx).successOutliers = groupSummaries(idx).meta.successOutliers;
    else
        logs(idx).successOutliers = struct('branchIdx',{},'attemptIndex',{},'value',{},'label',{},'primaryMetric',{},'secondaryMetric',{},'median',{},'scale',{},'deviation',{},'normalizedDeviation',{});
    end
end
end

function plotGroupedMetricFigure(cfg, lineSets, groupLabels)
if isempty(lineSets)
    return
end

figure(cfg.figureId); clf; hold on; grid on;
if isfield(cfg,'xlabel'), xlabel(cfg.xlabel); end
if isfield(cfg,'ylabel'), ylabel(cfg.ylabel); end
if isfield(cfg,'title'), title(cfg.title); end

uniqueGroups = unique(groupLabels, 'stable');
colors = lines(max(numel(uniqueGroups), 1));
lineStyles = {'-','--',':','-.'};
legendEntries = {};
legendHandles = [];
for idx = 1:numel(lineSets)
    line = lineSets(idx);
    if isempty(line.x)
        continue
    end
    groupIdx = find(strcmp(uniqueGroups, line.groupLabel), 1);
    if isempty(groupIdx)
        groupIdx = 1;
    end
    baseColor = colors(mod(groupIdx-1, size(colors,1))+1, :);
    color = branchColorVariant(baseColor, line.branchIdx);
    ls = lineStyles{mod(line.branchIdx-1, numel(lineStyles))+1};
    h = plot(line.x, line.y, 'LineWidth', 1.5, 'Color', color, 'LineStyle', ls);
    legendEntries{end+1} = line.lineLabel; %#ok<AGROW>
    legendHandles(end+1) = h; %#ok<AGROW>
end
 [legendHandles, legendEntries] = appendAsymptoteLegendEntries(cfg, lineSets, uniqueGroups, colors, legendHandles, legendEntries);
if ~isempty(legendEntries)
    legendEntries = cellfun(@formatLegendLabel, legendEntries, 'UniformOutput', false);
    if numel(legendHandles) == numel(legendEntries)
        legend(legendHandles, legendEntries,'Location','best'); box on;
    else
        legend(legendEntries,'Location','best'); box on;
    end
end
plotGroupedAsymptoteLines(cfg, lineSets);
applyNumericAxisFormatting(gca);
end

function value = getLineField(line, fieldName, defaultValue)
if isfield(line, fieldName) && ~isempty(line.(fieldName))
    value = line.(fieldName);
else
    value = defaultValue;
end
end

function plotGroupedAsymptoteLines(cfg, lineSets)
if isempty(lineSets)
    return
end
if ~shouldPlotGroupedAsymptote(cfg)
    return
end
ax = gca;
yl = get(ax, 'YLim');
if isempty(yl) || numel(yl) ~= 2
    return
end
groupLabels = {lineSets.groupLabel};
uniqueGroups = unique(groupLabels, 'stable');
colors = lines(max(numel(uniqueGroups), 1));
for idx = 1:numel(lineSets)
    if ~isfield(lineSets(idx),'asymptoteX') || isempty(lineSets(idx).asymptoteX)
        continue
    end
    groupIdx = find(strcmp(uniqueGroups, lineSets(idx).groupLabel), 1);
    if isempty(groupIdx)
        groupIdx = 1;
    end
    baseColor = colors(mod(groupIdx-1, size(colors,1))+1, :);
    color = branchColorVariant(baseColor, lineSets(idx).branchIdx);
    uniqX = unique(lineSets(idx).asymptoteX(isfinite(lineSets(idx).asymptoteX)));
    for k = 1:numel(uniqX)
        plot([uniqX(k) uniqX(k)], yl, '-.', 'LineWidth', 1.0, 'Color', color, 'HandleVisibility', 'off');
    end
end
end

% Objective: Add asymptote legend entries per M group.
% Purpose: Summarize discontinuities with consistent styling.
% SWOT: S-clear legend; W-more legend rows; O-reuse in grouped plots; T-clutter if many groups.
function [handlesOut, entriesOut] = appendAsymptoteLegendEntries(cfg, lineSets, uniqueGroups, colors, legendHandles, legendEntries)
handlesOut = legendHandles;
entriesOut = legendEntries;
if isempty(lineSets)
    return
end
axisLabel = sanitizeAxisLabel(getfieldWithDefault(cfg,'xlabel','x'));
groupLabels = {lineSets.groupLabel};
for g = 1:numel(uniqueGroups)
    groupLabel = uniqueGroups{g};
    idxList = find(strcmp(groupLabels, groupLabel));
    if isempty(idxList)
        continue
    end
    asymX = [];
    pickIdx = [];
    for k = 1:numel(idxList)
        line = lineSets(idxList(k));
        if isfield(line,'asymptoteX') && ~isempty(line.asymptoteX)
            asymX = [asymX; line.asymptoteX(:)]; %#ok<AGROW>
            if isempty(pickIdx)
                pickIdx = idxList(k);
            end
        end
    end
    if isempty(asymX)
        continue
    end
    asymX = unique(asymX(isfinite(asymX)));
    if isempty(asymX)
        continue
    end
    xsTxt = strjoin(arrayfun(@formatNumericToken, asymX, 'UniformOutput', false), ', ');
    label = sprintf('%s asymptote @ %s=%s', groupLabel, axisLabel, xsTxt);
    baseColor = colors(mod(g-1, size(colors,1))+1, :);
    branchIdx = getLineField(lineSets(pickIdx),'branchIdx',1);
    color = branchColorVariant(baseColor, branchIdx);
    h = plot(nan, nan, '-.', 'LineWidth', 1.0, 'Color', color, 'HandleVisibility', 'off');
    handlesOut(end+1) = h; %#ok<AGROW>
    entriesOut{end+1} = label; %#ok<AGROW>
end
end

function tf = shouldPlotGroupedAsymptote(cfg)
tf = false;
if nargin < 1 || ~isstruct(cfg)
    return
end
if ~isfield(cfg,'asymptote') || ~isstruct(cfg.asymptote) || ~isfield(cfg.asymptote,'enabled')
    return
end
tf = logical(cfg.asymptote.enabled);
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

function writeAttemptLogSection(fid, attemptLog, branchCount, label, figCfg)
if isempty(attemptLog)
    return
end
if nargin < 4 || strlength(string(label)) == 0
    sectionKey = 'attempt_log';
else
    sectionKey = sprintf('attempt_log_%s', sanitizeSectionLabel(label));
end
writeCsvRow(fid, {'section', sectionKey});

branchCount = max(branchCount, 1);
header = {'attempt_index','sweep_value','sweep_label','all_branches_success','branches_deviation','branches_deviation_cf','branches_deviation_nu'};
for b = 1:branchCount
    primaryLabel = buildMetricColumnLabel(attemptLog, b, 'primaryMetricLabels', 'primary_metric');
    secondaryLabel = buildMetricColumnLabel(attemptLog, b, 'secondaryMetricLabels', 'secondary_metric');
    header = [header, ...
        sprintf('B%d_status', b), ...
        sprintf('B%d_configured_step', b), ...
        sprintf('B%d_sweep_step', b), ...
        sprintf('B%d_avg_mesh_step', b), ...
        sprintf('B%d_max_residual', b), ...
        sprintf('B%d_error_message', b), ...
        sprintf('B%d_state', b), ...
        sprintf('B%d_fpp0', b), ...
        sprintf('B%d_neg_thp0', b), ...
        sprintf('B%d_initial_guess_error', b), ...
        sprintf('B%d_iterations', b), ...
        sprintf('B%d_%s', b, primaryLabel), ...
        sprintf('B%d_%s', b, secondaryLabel), ...
        sprintf('B%d_initial_guess_mesh', b), ...
        sprintf('B%d_initial_guess_profile', b)]; %#ok<AGROW>
end
writeCsvRow(fid, header);

for idx = 1:numel(attemptLog)
    attempt = attemptLog(idx);
    row = {idx, attempt.value, attempt.label, logical(getfieldWithDefault(attempt, 'allBranchesSuccess', false)), ...
        getfieldWithDefault(attempt, 'branchesDeviation', NaN), ...
        getfieldWithDefault(attempt, 'branchesDeviationPrimary', NaN), ...
        getfieldWithDefault(attempt, 'branchesDeviationSecondary', NaN)};
    statuses = attempt.branchStatus;
    cfgSteps = attempt.stepSizes;
    sweepSteps = attempt.sweepSteps;
    avgSteps = attempt.avgMeshSteps;
    residuals = attempt.maxResiduals;
    errors = attempt.errorMessages;
    guessErrors = attempt.initialGuessErrors;
    iterCounts = attempt.iterations;
    solVals = attempt.solutionValues;
    secVals = attempt.secondaryValues;
    meshCells = attempt.initialGuessMesh;
    profileCells = attempt.initialGuessProfiles;
    for b = 1:branchCount
        statusVal = '';
        if b <= numel(statuses) && ~isempty(statuses{b})
            statusVal = statuses{b};
        end
        cfgVal = pickAttemptValue(cfgSteps, b);
        sweepStep = pickAttemptValue(sweepSteps, b);
        avgVal = pickAttemptValue(avgSteps, b);
        resVal = pickAttemptValue(residuals, b);
        errVal = '';
        if b <= numel(errors)
            errVal = errors{b};
        end
        stateTxt = formatStateVectorAttempt(attempt, b);
        fpp0Val = extractStateComponent(attempt, b, 3, 1);
        negThp0Val = extractStateComponent(attempt, b, 5, -1);
        guessErr = pickAttemptValue(guessErrors, b);
        iterVal = pickAttemptValue(iterCounts, b);
        solVal = pickAttemptValue(solVals, b);
        secVal = pickAttemptValue(secVals, b);
        meshTxt = formatGuessArrayForExport(fetchAttemptCell(meshCells, b));
        profileTxt = formatGuessArrayForExport(fetchAttemptCell(profileCells, b));
        row = [row, {statusVal, cfgVal, sweepStep, avgVal, resVal, errVal, stateTxt, fpp0Val, negThp0Val, guessErr, iterVal, solVal, secVal, meshTxt, profileTxt}]; %#ok<AGROW>
    end
    writeCsvRow(fid, row);
end
writeCsvRow(fid, {''});
logAttemptEntries(figCfg, label, attemptLog, branchCount);
end

function writeContinuousRangeSection(fid, ranges, label, statusKey)
if isempty(ranges)
    return
end
if nargin < 3
    label = '';
end
if nargin < 4 || strlength(string(statusKey)) == 0
    statusKey = 'status';
end
statusKey = lower(strtrim(string(statusKey)));
if strlength(string(label)) == 0
    sectionKey = sprintf('continuous_%s_ranges', statusKey);
else
    sectionKey = sprintf('continuous_%s_ranges_%s', statusKey, sanitizeSectionLabel(label));
end
writeCsvRow(fid, {'section', sectionKey});
header = {'branch_index','start_attempt','end_attempt','count','start_value','end_value','start_label','end_label'};
writeCsvRow(fid, header);
for idx = 1:numel(ranges)
    entry = ranges(idx);
    startLabel = formatRangeLabel(entry.startLabel);
    endLabel = formatRangeLabel(entry.endLabel);
    row = {entry.branchIdx, entry.startAttempt, entry.endAttempt, entry.count, ...
        entry.startValue, entry.endValue, startLabel, endLabel};
    writeCsvRow(fid, row);
end
writeCsvRow(fid, {''});
end

function writeSuccessOutlierSection(fid, outliers, label)
if isempty(outliers)
    return
end
if nargin < 2 || isempty(outliers)
    return
end
if nargin < 3
    label = '';
end
if strlength(string(label)) == 0
    sectionKey = 'success_outliers';
else
    sectionKey = sprintf('success_outliers_%s', sanitizeSectionLabel(label));
end
writeCsvRow(fid, {'section', sectionKey});
header = {'attempt_index','branch','sweep_value','sweep_label','primary_metric','secondary_metric','deviation','normalized_deviation','median','scale'};
writeCsvRow(fid, header);
for idx = 1:numel(outliers)
    entry = outliers(idx);
    sweepLabel = formatRangeLabel(entry.label);
    row = {entry.attemptIndex, entry.branchIdx, entry.value, sweepLabel, entry.primaryMetric, entry.secondaryMetric, entry.deviation, entry.normalizedDeviation, entry.median, entry.scale};
    writeCsvRow(fid, row);
end
writeCsvRow(fid, {''});
end

function logAttemptEntries(figCfg, label, attemptLog, branchCount)
if isempty(attemptLog) || branchCount < 1
    return
end
if nargin < 1 || isempty(figCfg)
    figName = 'figure';
else
    if isfield(figCfg,'title') && ~isempty(figCfg.title)
        figName = figCfg.title;
    elseif isfield(figCfg,'figureId')
        figName = sprintf('Figure %d', figCfg.figureId);
    else
        figName = 'figure';
    end
end
if nargin < 2 || strlength(string(label)) == 0
    groupSuffix = '';
else
    groupSuffix = sprintf(' | %s', label);
end
xAxisLabel = sanitizeAxisLabel(getfieldWithDefault(figCfg,'xlabel','x'));
for idx = 1:numel(attemptLog)
    attempt = attemptLog(idx);
    sweepVal = attempt.value;
    sweepLabel = attempt.label;
    % for b = 1:branchCount
    %     lineLogIdx = fetchNextLineLogIndex();
    %     emitFactoryLog('[L%dS] -----------------------\n', lineLogIdx);
    %     emitFactoryLog('=======================\n');
    %     statusVal = '';
    %     if b <= numel(attempt.branchStatus) && ~isempty(attempt.branchStatus{b})
    %         statusVal = attempt.branchStatus{b};
    %     end
    %     primaryVal = pickAttemptValue(attempt.solutionValues, b);
    %     secondaryVal = pickAttemptValue(attempt.secondaryValues, b);
    %     cfgStep = pickAttemptValue(attempt.stepSizes, b);
    %     sweepStep = pickAttemptValue(attempt.sweepSteps, b);
    %     consoleTxt = '';
    %     if b <= numel(attempt.consoleLogs) && ~isempty(attempt.consoleLogs{b})
    %         consoleTxt = attempt.consoleLogs{b};
    %     end
    %     emitFactoryLog('Figure "%s"%s | branch %d | sweep=%g (%s) [delta=%.6g] | %s-configured_step=%g | primary=%.6g | secondary=%.6g | status=%s\n', ...
    %         figName, groupSuffix, b, sweepVal, formatSweepLabel(sweepLabel), sweepStep, xAxisLabel, cfgStep, primaryVal, secondaryVal, statusVal);
    %     if ~isempty(strtrim(consoleTxt))
    %         emitFactoryLog('%s\n', strtrim(consoleTxt));
    %     end
    %     emitFactoryLog('=======================\n');
    %     emitFactoryLog('[L%dE] -----------------------\n', lineLogIdx);
    % end
end
end

function txt = formatSweepLabel(label)
if isstring(label)
    txt = char(label);
elseif ischar(label)
    txt = label;
else
    txt = char(string(label));
end
if isempty(txt)
    txt = 'n/a';
end
end

function txt = formatRangeLabel(label)
if isstring(label)
    txt = char(label);
elseif ischar(label)
    txt = label;
else
    txt = char(string(label));
end
end

function token = sanitizeSectionLabel(label)
txt = lower(strtrim(string(label)));
if strlength(txt) == 0
    token = 'section';
    return
end
txt = regexprep(txt, '[^a-z0-9]+', '_');
txt = regexprep(txt, '_+', '_');
txt = regexprep(txt, '^_|_$', '');
if strlength(txt) == 0
    token = 'section';
else
    token = char(txt);
end
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

function value = fetchAttemptCell(values, idx)
value = [];
if idx < 1
    return
end
if iscell(values)
    if idx <= numel(values)
        value = values{idx};
    end
elseif ~isempty(values)
    if idx <= numel(values)
        value = values(idx);
    end
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
function label = buildMetricColumnLabel(attemptLog, branchIdx, fieldName, defaultLabel)
if nargin < 4 || isempty(defaultLabel)
    defaultLabel = 'metric';
end
rawLabel = extractMetricLabel(attemptLog, branchIdx, fieldName);
if strlength(string(rawLabel)) == 0
    rawLabel = defaultLabel;
end
label = sanitizeMetricLabel(rawLabel);
end

function label = extractMetricLabel(attemptLog, branchIdx, fieldName)
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

function label = sanitizeMetricLabel(labelIn)
text = lower(strtrim(string(labelIn)));
if strlength(text) == 0
    label = 'metric';
    return
end
text = regexprep(text, '[^a-z0-9]+', '_');
text = regexprep(text, '_+', '_');
text = regexprep(text, '^_|_$', '');
if strlength(text) == 0
    label = 'metric';
else
    label = char(text);
end
end

function label = deriveOutputLabel(cfg, sweepSummaries)
if isfield(cfg,'outputTag') && strlength(string(cfg.outputTag)) > 0
    label = sanitizeRunToken(cfg.outputTag);
    return
end
for idx = 1:numel(sweepSummaries)
    entry = sweepSummaries(idx);
    if ~isstruct(entry) || ~isfield(entry,'meta') || isempty(entry.meta)
        continue
    end
    meta = entry.meta;
    if ~isfield(meta,'sweepValues') || isempty(meta.sweepValues)
        continue
    end
    finiteVals = meta.sweepValues(isfinite(meta.sweepValues));
    if isempty(finiteVals)
        continue
    end
    sweepName = meta.sweepName;
    if strlength(string(sweepName)) == 0
        sweepName = sprintf('sweep%d', idx);
    end
    label = sprintf('%s_%s', sanitizeRunToken(sweepName), buildRangeLabel(finiteVals));
    label = sanitizeRunToken(label);
    return
end
label = 'run';
end

function token = sanitizeRunToken(value)
strVal = strtrim(char(string(value)));
if isempty(strVal)
    token = 'run';
    return
end
strVal = lower(strVal);
strVal = regexprep(strVal, '[^a-z0-9]+', '_');
strVal = regexprep(strVal, '_+', '_');
strVal = regexprep(strVal, '^_|_$', '');
if isempty(strVal)
    token = 'run';
else
    token = strVal;
end
end

function window = resolveFigureWindow(cfg)
window = [];
if isfield(cfg,'domainWindow') && ~isempty(cfg.domainWindow)
    window = cfg.domainWindow;
end
end

function label = sanitizeAxisLabel(value)
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

function val = getfieldWithDefault(s, fieldName, defaultVal)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    val = s.(fieldName);
else
    val = defaultVal;
end
end

function cfg = ensureDomainDefaults(cfg)
if ~isfield(cfg,'domainMinList') || isempty(cfg.domainMinList)
    cfg.domainMinList = 0;
end
if ~isfield(cfg,'domainMaxList') || isempty(cfg.domainMaxList)
    cfg.domainMaxList = 1;
end
if ~isfield(cfg,'domainGridSizeList') || isempty(cfg.domainGridSizeList)
    cfg.domainGridSizeList = 401;
end
if ~isfield(cfg,'domainLabel') || isempty(cfg.domainLabel)
    cfg.domainLabel = 'x';
end
if ~isfield(cfg,'numericFormat') || ~isstruct(cfg.numericFormat)
    cfg.numericFormat = struct();
end
if ~isfield(cfg.numericFormat,'decimals') || isempty(cfg.numericFormat.decimals)
    cfg.numericFormat.decimals = 8;
end
if ~isfield(cfg.numericFormat,'removeLeadingZero') || isempty(cfg.numericFormat.removeLeadingZero)
    cfg.numericFormat.removeLeadingZero = true;
end
if ~isfield(cfg.numericFormat,'removeTrailingZeros') || isempty(cfg.numericFormat.removeTrailingZeros)
    cfg.numericFormat.removeTrailingZeros = false;
end
end

function cfg = applyMethodDomainConfig(cfg, method)
if ~isstruct(cfg)
    cfg = struct();
end
if nargin < 2 || ~isstruct(method)
    return
end
if isfield(method,'domainGridSizeList') && ~isempty(method.domainGridSizeList)
    if ~isfield(cfg,'domainGridSizeList') || isempty(cfg.domainGridSizeList)
        cfg.domainGridSizeList = method.domainGridSizeList;
    end
end
end

function method = ensureMethodDefaults(method, cfg)
if nargin < 1 || isempty(method)
    method = struct();
end
if nargin < 2
    cfg = struct();
end
if ~isfield(method,'solver') || strlength(string(method.solver)) == 0
    method.solver = 'bvp4c';
end
if ~isfield(method,'displayName') || strlength(string(method.displayName)) == 0
    method.displayName = upper(method.solver);
end
method.relativeTolerance = getfieldWithDefault(method,'relativeTolerance',1e-6);
method.absoluteTolerance = getfieldWithDefault(method,'absoluteTolerance',1e-8);
method.maxMeshPoints = getfieldWithDefault(method,'maxMeshPoints',60000);
method.maxSolverAttempts = getfieldWithDefault(method,'maxSolverAttempts',getfieldWithDefault(cfg,'maxSolverAttempts',10));
method.timeLimitSeconds = getfieldWithDefault(method,'timeLimitSeconds',getfieldWithDefault(cfg,'timeLimitSeconds',10));
method.attemptTimeLimit = getfieldWithDefault(method,'attemptTimeLimit',min(10, method.timeLimitSeconds));
method.disableContinuation = getfieldWithDefault(method,'disableContinuation',getfieldWithDefault(cfg,'disableContinuation',false));
domainLength = resolveDomainLength(cfg);
defaultStep = domainLength / 200;
method.initialStepSize = getfieldWithDefault(method,'initialStepSize',defaultStep);
method.minStepSize = getfieldWithDefault(method,'minStepSize',defaultStep/1000);
method.maxStepSize = getfieldWithDefault(method,'maxStepSize',domainLength/10);
method.stepSafetyFactor = getfieldWithDefault(method,'stepSafetyFactor',0.9);
method.stepGrowthLimit = getfieldWithDefault(method,'stepGrowthLimit',4);
method.stepShrinkLimit = getfieldWithDefault(method,'stepShrinkLimit',0.25);
method.maxStepCount = getfieldWithDefault(method,'maxStepCount',20000);
method.newtonTolerance = getfieldWithDefault(method,'newtonTolerance',1e-8);
method.newtonMaxIterations = getfieldWithDefault(method,'newtonMaxIterations',12);
method.jacobianPerturbation = getfieldWithDefault(method,'jacobianPerturbation',1e-6);
method.shootingTolerance = getfieldWithDefault(method,'shootingTolerance',1e-7);
method.shootingMaxIterations = getfieldWithDefault(method,'shootingMaxIterations',10);
method.shootingPerturbation = getfieldWithDefault(method,'shootingPerturbation',1e-6);
method.freeInitialIndices = getfieldWithDefault(method,'freeInitialIndices',[]);
method.initialStateProjector = getfieldWithDefault(method,'initialStateProjector',[]);
method.activeResidualIndices = getfieldWithDefault(method,'activeResidualIndices',[]);
method.integratorMaxRelTol = getfieldWithDefault(method,'integratorMaxRelTol',1e-4);
method.integratorMaxAbsTol = getfieldWithDefault(method,'integratorMaxAbsTol',1e-6);
method.integratorRelaxFactor = getfieldWithDefault(method,'integratorRelaxFactor',5);
end

function len = resolveDomainLength(cfg)
len = 1;
if nargin == 0 || isempty(cfg)
    return
end
minVal = pickValue(cfg.domainMinList, 1);
maxVal = pickValue(cfg.domainMaxList, 1);
if isnumeric(minVal) && isnumeric(maxVal)
    span = maxVal - minVal;
    if isfinite(span) && span > 0
        len = span;
    end
end
end

function stepSize = computePointSteps(xVals)
n = numel(xVals);
stepSize = nan(n,1);
if n <= 1
    return
end
stepSize(2:end) = diff(xVals);
end

function writeCsvRow(fid, cells)
if isempty(cells) || (numel(cells) == 1 && (isstring(cells{1}) || ischar(cells{1})) && strlength(string(cells{1})) == 0)
    fprintf(fid, '\n');
    return
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

function txt = formatStateVectorAttempt(attempt, branchIdx)
stateTxt = 'N/A';
if isfield(attempt,'stateVectors') && ~isempty(attempt.stateVectors) && branchIdx <= numel(attempt.stateVectors)
    vec = attempt.stateVectors{branchIdx};
    if ~isempty(vec)
        if isstring(vec) || ischar(vec)
            stateTxt = char(vec);
        elseif isnumeric(vec)
            stateTxt = formatNumericArray(vec);
        end
    end
end
txt = stateTxt;
end

function token = formatNumericToken(value)
token = numericFormat('token', value);
end

function txt = formatNumericArray(values)
txt = numericFormat('array', values);
end

function applyNumericFormatOptions(cfg)
numericFormat('apply', cfg);
end

function val = extractStateComponent(attempt, branchIdx, compIdx, scale)
val = NaN;
if nargin < 4 || isempty(scale)
    scale = 1;
end
if ~isfield(attempt,'stateVectors') || isempty(attempt.stateVectors)
    return
end
if branchIdx > numel(attempt.stateVectors)
    return
end
vec = attempt.stateVectors{branchIdx};
if isempty(vec) || compIdx > numel(vec)
    return
end
val = scale * vec(compIdx);
end

function emitFactoryLog(fmt, varargin)
txt = sprintf(fmt, varargin{:});
fprintf('%s', txt);
if isappdata(0,'pehf_logPath')
    logPath = getappdata(0,'pehf_logPath');
    if ~isempty(logPath)
        try
            fid = fopen(logPath,'a');
            if fid >= 0
                fprintf(fid,'%s', txt);
                fclose(fid);
            end
        catch
        end
    end
end
end

function idx = fetchNextFigureLogIndex()
idx = fetchAppCounter('pehf_figureLogCounter');
end

function idx = fetchNextLineLogIndex()
idx = fetchAppCounter('pehf_lineLogCounter');
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

function cfg = updateSeedLibrary(cfg, sweepSummary)
if nargin < 2 || isempty(sweepSummary)
    return
end
if ~isfield(cfg,'seedLibrary')
    cfg.seedLibrary = struct();
end
if ~isfield(sweepSummary,'meta') || ~isfield(sweepSummary.meta,'sweepValues')
    return
end
values = sweepSummary.meta.sweepValues;
if isempty(values)
    return
end
sweepName = '';
if isfield(sweepSummary,'name')
    sweepName = sweepSummary.name;
end
if isfield(sweepSummary.meta,'sweepName') && ~isempty(sweepSummary.meta.sweepName)
    sweepName = sweepSummary.meta.sweepName;
end
key = sanitizeSeedKey(sweepName);
if strlength(key) == 0
    return
end
results = sweepSummary.results;
if isempty(results)
    return
end
branchCount = size(results,1);
if ~isfield(cfg.seedLibrary, key) || numel(cfg.seedLibrary.(key)) ~= branchCount
    cfg.seedLibrary.(key) = repmat(struct('entries',[]), branchCount, 1);
end
for b = 1:branchCount
    for k = 1:numel(values)
        sol = results(b,k).sol;
        if isempty(sol)
            continue
        end
        entry = struct('value', values(k), 'sol', sol);
        cfg.seedLibrary.(key)(b).entries = [cfg.seedLibrary.(key)(b).entries, entry]; %#ok<AGROW>
    end
end
end

function cfgOut = applySeedGuesses(cfgRoot, cfgOut, paramName, paramValue)
cfgOut = ensureDomainDefaults(cfgOut);
if nargin < 3 || strlength(string(paramName)) == 0
    return
end
if isfield(cfgOut,'disableSeedGuesses') && cfgOut.disableSeedGuesses
    return
end
if ~isfield(cfgRoot,'seedLibrary')
    return
end
key = sanitizeSeedKey(paramName);
if strlength(key) == 0 || ~isfield(cfgRoot.seedLibrary, key)
    return
end
entriesPerBranch = cfgRoot.seedLibrary.(key);
if isempty(entriesPerBranch)
    return
end
if ~iscell(cfgOut.guesses) || isempty(cfgOut.guesses)
    cfgOut.guesses = cfgRoot.guesses;
end
tol = max(1e-6, 1e-6 * max(1, abs(paramValue)));
for b = 1:min(numel(entriesPerBranch), numel(cfgOut.guesses))
    entries = entriesPerBranch(b).entries;
    if isempty(entries)
        continue
    end
    values = [entries.value];
    diffs = abs(values - paramValue);
    [minDiff, idxMin] = min(diffs);
    if isempty(idxMin) || minDiff > tol
        continue
    end
    sol = entries(idxMin).sol;
    cfgOut.guesses{b} = @(coord) evalSolution(sol, coord);
end
end

function key = sanitizeSeedKey(name)
txt = lower(strtrim(string(name)));
if strlength(txt) == 0
    key = '';
    return
end
txt = regexprep(txt, '[^a-z0-9]+', '');
key = char(txt);
end

function out = quoteString(str)
if isempty(str)
    out = '""';
else
    out = ['"', strrep(str, '"', '""'), '"'];
end
end

function val = pickValue(list, idx)
if nargin < 2 || isempty(idx)
    idx = 1;
end
if isempty(list)
    val = NaN;
    return
end
if isscalar(list)
    val = list;
    return
end
idx = max(1, min(numel(list), idx));
val = list(idx);
end

function yq = evalSolution(sol, xq)
persistent canUseDeval
if isempty(canUseDeval)
    canUseDeval = true;
end

if canUseDeval
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

xq = xq(:).';
yInterp = interp1(sol.x, sol.y.', xq, 'pchip', 'extrap').';

if isnumeric(xq) && isscalar(xq)
    yq = yInterp(:);
else
    yq = yInterp;
end
end
