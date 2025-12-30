function ranges = computeContinuousFailRanges(attemptLog, targetStatus)
%COMPUTECONTINUOUSFAILRANGES Identify consecutive sequences per branch.
%
% targetStatus (optional) : status string to match (default 'fail')
%
% ranges: struct array with fields
%   branchIdx, startAttempt, endAttempt, count,
%   startValue, endValue, startLabel, endLabel

ranges = struct('branchIdx',{},'startAttempt',{},'endAttempt',{}, ...
    'count',{},'startValue',{},'endValue',{},'startLabel',{},'endLabel',{});

if nargin == 0 || isempty(attemptLog)
    return
end
if nargin < 2 || strlength(string(targetStatus)) == 0
    targetStatus = 'fail';
end
targetStatus = lower(strtrim(string(targetStatus)));

branchCount = 0;
for idx = 1:numel(attemptLog)
    branchCount = max(branchCount, numel(attemptLog(idx).branchStatus));
end
if branchCount == 0
    return
end

for branchIdx = 1:branchCount
    active = struct('startIdx',[], 'startValue', NaN, 'startLabel', "", 'endIdx', [], 'endValue', NaN, 'endLabel', "");
    for attemptIdx = 1:numel(attemptLog)
        statusVal = '';
        entry = attemptLog(attemptIdx);
        if branchIdx <= numel(entry.branchStatus) && ~isempty(entry.branchStatus{branchIdx})
            statusVal = lower(strtrim(string(entry.branchStatus{branchIdx})));
        end
        isMatch = strcmp(statusVal, targetStatus);
        if isMatch
            if isempty(active.startIdx)
                active.startIdx = attemptIdx;
                active.startValue = fetchScalar(entry.value);
                active.startLabel = string(entry.label);
            end
            active.endIdx = attemptIdx;
            active.endValue = fetchScalar(entry.value);
            active.endLabel = string(entry.label);
        else
            if ~isempty(active.startIdx)
                ranges(end+1) = finalizeRangeStruct(branchIdx, active); %#ok<AGROW>
                active = struct('startIdx',[], 'startValue', NaN, 'startLabel', "", 'endIdx', [], 'endValue', NaN, 'endLabel', "");
            end
        end
    end
    if ~isempty(active.startIdx)
        ranges(end+1) = finalizeRangeStruct(branchIdx, active); %#ok<AGROW>
    end
end
end

function val = fetchScalar(value)
if isnumeric(value) && isscalar(value)
    val = value;
else
    val = NaN;
end
end

function entry = finalizeRangeStruct(branchIdx, active)
entry = struct( ...
    'branchIdx', branchIdx, ...
    'startAttempt', active.startIdx, ...
    'endAttempt', active.endIdx, ...
    'count', active.endIdx - active.startIdx + 1, ...
    'startValue', active.startValue, ...
    'endValue', active.endValue, ...
    'startLabel', active.startLabel, ...
    'endLabel', active.endLabel);
end
