function cfg = config()
% Objective: Centralize plotting/output configuration for the Arif profile.
% Purpose: Keep render and formatting options out of values.m.
% SWOT: S-clear separation; W-more indirection; O-reuse by other profiles; T-drift from values.m defaults.

cfg = struct();
cfg.plotOutliers = false;
cfg.disableSeedGuesses = false;
cfg.maxSolverAttempts = 2;
cfg.timeLimitSeconds = 2;
cfg.disableContinuation = true;

cfg.numericFormat = struct( ...
    'decimals', 8, ...
    'removeLeadingZero', false, ...
    'removeTrailingZeros', true);
end
