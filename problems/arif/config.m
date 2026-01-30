function cfg = config()
% Objective: Centralize plotting/output configuration for the Arif profile.
% Purpose: Keep render and formatting options out of values.m.
% SWOT: S-clear separation; W-more indirection; O-reuse by other profiles; T-drift from values.m defaults.

cfg = struct();
cfg.plotOutliers = true;
% Outlier filtering: "global" hides any point flagged on any branch; "branch" hides only per-branch outliers.
cfg.outlierFilterMode = 'global';
% Outlier sigma threshold: points beyond sigma * (robust spread) from median are flagged (higher = fewer outliers).
cfg.outlierSigma = 4.5;
% Toggle plotting points where all branches have zero deviation (true = keep, false = hide).
cfg.plotZeroDeviationPoints = true;
% Absolute tolerance for treating branch deviation as zero.
cfg.zeroDeviationTolerance = 1e-6;
cfg.disableSeedGuesses = false;
cfg.maxSolverAttempts = 2;
cfg.timeLimitSeconds = 4;
cfg.disableContinuation = true;
% Enable smart lambda spacing (uses values.sweep.smartLinspaceLambdaVals).
% cfg.smartLinspace = true;

cfg.numericFormat = struct( ...
    'decimals', 8, ...
    'removeLeadingZero', false, ...
    'removeTrailingZeros', true);
end
