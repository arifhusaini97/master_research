% ========================Start BVP4C==========================================
% function spec = method(cfg)
% % METHOD Numerical strategy for Arif profile (bvp4c collocation).
%
%     if nargin < 1
%         cfg = struct();
%     end
%
%     spec = struct();
%     spec.solver = 'bvp4c';
%     spec.displayName = 'Collocation (bvp4c)';
%     spec.relativeTolerance = 1e-10;
%     spec.absoluteTolerance = 1e-10;
%     spec.maxMeshPoints = 60000;
%     spec.maxSolverAttempts = 10;
%     spec.timeLimitSeconds = 10;
% end
% ========================End BVP4C==========================================
% ========================Start BVP5C==========================================
function spec = method(cfg)
% METHOD Numerical strategy for Arif profile (bvp4c collocation).

if nargin < 1
    cfg = struct();
end

spec = struct();
spec.solver = 'bvp5c';
spec.displayName = 'Collocation (bvp5c)';
spec.relativeTolerance = 1e-10;
spec.absoluteTolerance = 1e-10;
spec.maxMeshPoints = 60000;
spec.maxSolverAttempts = 5;
spec.timeLimitSeconds = 5;
spec.meshLabel = '\eta';
end
% ========================End BVP5C==========================================
% ========================Start BDF2==========================================
% function spec = method(cfg)
% % METHOD Numerical strategy definition for Natasha (BDF2-first workflow).
% % Objective: Configure the shooting-based BDF2 engine for Natasha''s stiff profile.
% % Purpose : Encapsulate per-problem tolerances so the engine can remain generic.
% % SWOT     S: tailored convergence knobs; W: assumes valid cfg.p; O: reusable template for other BDF cases; T: projector must stay in sync with BCs.
%
%     if nargin < 1
%         cfg = struct();
%     end
%     if ~isfield(cfg, 'p')
%         cfg.p = struct('S', 1, 'lambda', 0, 'L1', 0, 'L2', 0);
%     end
%
%     spec = struct();
%     spec.solver = 'bdf2';
%     spec.displayName = 'Adaptive BDF2 shooter';
%     spec.relativeTolerance = 1e-8;
%     spec.absoluteTolerance = 1e-10;
%     spec.initialStepSize = 1e-2;
%     spec.minStepSize = 1e-4;
%     spec.maxStepSize = 0.12;
%     spec.stepSafetyFactor = 0.8;
%     spec.stepGrowthLimit = 1.6;
%     spec.stepShrinkLimit = 0.35;
%     spec.maxStepCount = 150000;
%     spec.frontierRefineLimit = 0.25;
%     spec.frontierRefineStep = 4e-4;
%     spec.maxSolverAttempts = 1;
%     spec.timeLimitSeconds = 10;
%     spec.attemptTimeLimit = 6;
%     spec.shootingTolerance = 1e-7;
%     spec.shootingMaxIterations = 8;
%     spec.shootingPerturbation = 5e-5;
%     spec.freeInitialIndices = [3, 5];
%     spec.activeResidualIndices = [4, 5];
%     spec.seedWithCollocation = true;
%     spec.seedCollocationRelTol = 1e-6;
%     spec.seedCollocationAbsTol = 1e-8;
%     spec.seedCollocationMaxPoints = 10000;
%
%     tSlopeMag = 35;
%     spec.initialThetaSlope = tSlopeMag;
%     spec.initialStateProjector = @(ya) projectNatashaInitialState(ya, cfg.p, tSlopeMag);
% end
%
% function yaProj = projectNatashaInitialState(ya, params, tSlopeMag)
%     yaProj = ya(:);
%     yaProj(1) = params.S;
%     yaProj(4) = 1;
%     baseArg = max(1 - params.L2 * yaProj(3), 1e-6);
%     yaProj(2) = params.lambda + params.L1 * yaProj(3) / sqrt(baseArg);
%     if yaProj(5) >= 0
%         yaProj(5) = -max(tSlopeMag, yaProj(5));
%     elseif abs(yaProj(5)) < 1e-3
%         yaProj(5) = -1e-3;
%     end
% end
% ========================End BDF2==========================================
