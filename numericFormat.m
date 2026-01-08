function varargout = numericFormat(action, varargin)
% Objective: Centralize numeric formatting used in logs/CSV/labels.
% Purpose: Single source of truth for formatting options and tokens.
% SWOT: S-consistent formatting; W-extra indirection; O-reuse across modules; T-misuse of action keys.

if nargin == 0 || strlength(string(action)) == 0
    error('numericFormat:missingAction','An action name is required.');
end

switch lower(strtrim(string(action)))
    case "apply"
        if nargin < 2
            error('numericFormat:missingConfig','apply requires a config struct.');
        end
        applyNumericFormatOptions(varargin{1});
        return
    case "get"
        varargout{1} = getNumericFormatOptions();
        return
    case "token"
        if nargin < 2
            error('numericFormat:missingValue','token requires a numeric value.');
        end
        varargout{1} = formatNumericToken(varargin{1});
        return
    case "array"
        if nargin < 2
            error('numericFormat:missingArray','array requires a numeric array.');
        end
        varargout{1} = formatNumericArray(varargin{1});
        return
    otherwise
        error('numericFormat:unknownAction','Unknown numericFormat action "%s".', action);
end
end

function applyNumericFormatOptions(cfg)
opts = getNumericFormatOptions();
if isstruct(cfg) && isfield(cfg,'numericFormat') && isstruct(cfg.numericFormat)
    opts = mergeNumericFormatOptions(opts, cfg.numericFormat);
end
setappdata(0, 'pehf_numericFormat', opts);
end

function token = formatNumericToken(value)
opts = getNumericFormatOptions();
decimals = opts.decimals;
scale = 10^decimals;
rounded = round(value * scale) / scale;
if abs(rounded) < (0.5 / scale)
    rounded = 0;
end
token = sprintf(['%.', num2str(decimals), 'f'], rounded);
if opts.removeTrailingZeros
    token = regexprep(token, '0+$', '');
    token = regexprep(token, '\.$', '');
end
if opts.removeLeadingZero
    token = regexprep(token, '^-0\.', '-.');
    token = regexprep(token, '^0\.', '.');
end
end

function txt = formatNumericArray(values)
vals = values(:).';
formatted = arrayfun(@formatNumericToken, vals, 'UniformOutput', false);
txt = ['[', strjoin(formatted, ', '), ']'];
end

function opts = getNumericFormatOptions()
opts = struct('decimals', 8, 'removeLeadingZero', true, 'removeTrailingZeros', false);
if isappdata(0, 'pehf_numericFormat')
    stored = getappdata(0, 'pehf_numericFormat');
    if isstruct(stored)
        opts = mergeNumericFormatOptions(opts, stored);
    end
end
opts = normalizeNumericFormatOptions(opts);
end

function opts = mergeNumericFormatOptions(base, update)
opts = base;
if isfield(update, 'decimals')
    opts.decimals = update.decimals;
end
if isfield(update, 'removeLeadingZero')
    opts.removeLeadingZero = update.removeLeadingZero;
end
if isfield(update, 'removeTrailingZeros')
    opts.removeTrailingZeros = update.removeTrailingZeros;
end
end

function opts = normalizeNumericFormatOptions(opts)
if ~isfield(opts, 'decimals') || isempty(opts.decimals) || ~isfinite(opts.decimals)
    opts.decimals = 8;
end
opts.decimals = max(0, round(opts.decimals));
if ~isfield(opts, 'removeLeadingZero') || isempty(opts.removeLeadingZero)
    opts.removeLeadingZero = true;
end
opts.removeLeadingZero = logical(opts.removeLeadingZero);
if ~isfield(opts, 'removeTrailingZeros') || isempty(opts.removeTrailingZeros)
    opts.removeTrailingZeros = false;
end
opts.removeTrailingZeros = logical(opts.removeTrailingZeros);
end
