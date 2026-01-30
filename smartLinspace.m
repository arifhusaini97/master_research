function cfg = smartLinspace(minVal, maxVal, count, minRate, maxRate)
%SMARTLINSPACE Describe slope-constrained sweep configuration.
% minRate/maxRate represent bounds on |delta(metric)/delta(lambda)|.

if nargin < 5
    error('smartLinspace:InvalidInput', 'Expected minVal, maxVal, count, minRate, maxRate.');
end
if ~isfinite(minVal) || ~isfinite(maxVal) || ~isfinite(count) ...
        || ~isfinite(minRate) || ~isfinite(maxRate)
    error('smartLinspace:InvalidInput', 'Inputs must be finite.');
end
count = round(count);
if count < 2
    count = 2;
end
if minRate < 0 || maxRate < 0 || maxRate < minRate
    error('smartLinspace:InvalidInput', 'Rate bounds must be non-negative and maxRate >= minRate.');
end

cfg = struct( ...
    'min', minVal, ...
    'max', maxVal, ...
    'count', count, ...
    'minRate', minRate, ...
    'maxRate', maxRate);
end
