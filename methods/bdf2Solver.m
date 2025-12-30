function [sol, diagnostics] = bdf2Solver(odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel, utils)
% Objective: Provide an adaptive BDF2 shooting integrator that mirrors the Celaya/Surati scheme.
% Purpose : Keeps each problem's tolerance/timeout knobs centralized while letting displayBaseFn stay agnostic.
% SWOT     S: reusable adaptive integrator; W: more complex than MATLAB's collocation; O: unlocks solver choice per problem; T: stiff IVPs may still need fine-tuning.

    diagnostics = utils.initBranchDiagnostics(stepSize);
    diagnostics.warningId = '';
    diagnostics.warningMsg = '';
    diagnostics.errorMessage = '';
    diagnostics.status = 'fail';
    diagnostics.attempts = 0;
    sol = [];

    startClock = tic;
    timeLimit = methodCfg.timeLimitSeconds;
    attemptLimitDefault = methodCfg.attemptTimeLimit;

    try
        ya = ensureColumnGuess(guessFcn, domainMin);
    catch err
        diagnostics.errorMessage = err.message;
        sol = [];
        return
    end

    if isfield(methodCfg,'seedWithCollocation') && methodCfg.seedWithCollocation
        ya = seedInitialStateWithCollocation(ya, odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel, utils);
    end

    if ~isempty(methodCfg.initialStateProjector) && isa(methodCfg.initialStateProjector,'function_handle')
        ya = methodCfg.initialStateProjector(ya);
    end
    stateDim = numel(ya);
    freeIdx = methodCfg.freeInitialIndices;
    if isempty(freeIdx)
        freeIdx = 1:stateDim;
    end

    shootTol = methodCfg.shootingTolerance;
    maxShootIter = methodCfg.shootingMaxIterations;
    maxShootIter = max(1, round(maxShootIter));
    attemptCap = utils.getfieldWithDefault(methodCfg,'maxSolverAttempts',maxShootIter);
    if isfinite(attemptCap) && attemptCap > 0
        maxShootIter = min(maxShootIter, round(attemptCap));
    end
    perturbation = methodCfg.shootingPerturbation;

    totalBcEvals = 0;
    bestResidual = Inf;
    bestSol = [];
    bestInfo = struct();
    stagnantCount = 0;

    for iter = 1:maxShootIter
        bcEvalCheckpoint = totalBcEvals;
        diagnostics.attempts = iter;
        utils.enforceTimeLimit(@(~)0, startClock, timeLimit); %#ok<NASGU>
        timeLeftGlobal = max(0, timeLimit - toc(startClock));
        attemptLimit = min(attemptLimitDefault, timeLeftGlobal);
        if attemptLimit <= 0
            diagnostics.errorMessage = 'Time limit exceeded for BDF2 solver.';
            sol = bestSol;
            return
        end
        attemptStart = tic;
        try
            [traj, info] = integrateBdf2Trajectory(odeFun, domainMin, domainMax, ya, methodCfg, attemptStart, attemptLimit, utils);
        catch err
            if strcmp(err.identifier,'displayBaseFn:timeout')
                diagnostics.errorMessage = err.message;
                diagnostics.warningId = err.identifier;
                diagnostics.warningMsg = err.message;
                sol = bestSol;
                return
            else
                rethrow(err);
            end
        end
        if info.success && ~isempty(traj.x)
            yb = traj.y(:,end);
            residualFull = bcFun(ya, yb);
            totalBcEvals = totalBcEvals + 1;
        else
            residualFull = Inf(stateDim,1);
        end
        if isempty(methodCfg.activeResidualIndices)
            residual = residualFull;
        else
            residual = residualFull(methodCfg.activeResidualIndices);
        end
        attemptBcEvals = max(0, totalBcEvals - bcEvalCheckpoint);
        info.bcEvals = attemptBcEvals;
        resNorm = max(abs(residual));
        if info.success && ~isempty(traj.x)
            logBody = utils.buildMethodConsoleLog(numel(traj.x), resNorm, info.fevals, attemptBcEvals);
            meta = methodCfg;
            meta.solutionRange = [traj.x(1), traj.x(end)];
            liveLog = utils.formatSolverConsoleLog(logBody, iter, domainMin, domainMax, domainLabel, meta);
            if ~isempty(liveLog)
                utils.emitRunLog('%s\n', liveLog);
            end
        else
            liveLog = '';
        end
        if resNorm < bestResidual
            bestResidual = resNorm;
            bestSol = traj;
            bestInfo = info;
            stagnantCount = 0;
        else
            stagnantCount = stagnantCount + 1;
        end
        if stagnantCount >= 2 && resNorm > shootTol
            break
        end
        if resNorm <= shootTol
            sol = traj;
            diagnostics.status = 'success';
            diagnostics.meshPoints = numel(traj.x);
            diagnostics.avgMeshStep = utils.averageMeshSpacing(traj.x);
            diagnostics.iterations = info.totalNewtonIters;
            diagnostics.initialGuessError = utils.computeInitialGuessError(sol, guessFcn);
            diagnostics.maxResidual = resNorm;
            if isempty(liveLog)
                logBody = utils.buildMethodConsoleLog(numel(traj.x), resNorm, info.fevals, info.bcEvals);
                meta = methodCfg;
                meta.solutionRange = [traj.x(1), traj.x(end)];
                diagnostics.consoleLog = utils.formatSolverConsoleLog(logBody, iter, domainMin, domainMax, domainLabel, meta);
            else
                diagnostics.consoleLog = liveLog;
            end
            return
        end

        nEq = numel(residual);
        J = zeros(nEq, numel(freeIdx));
        for col = 1:numel(freeIdx)
            idxVar = freeIdx(col);
            baseDelta = perturbation * max(1, abs(ya(idxVar)));
            delta = baseDelta;
            minDelta = max(baseDelta * 1e-6, eps);
            successCol = false;
            rPert = [];
            while delta >= minDelta
                yaPert = ya;
                yaPert(idxVar) = yaPert(idxVar) + delta;
                if ~isempty(methodCfg.initialStateProjector) && isa(methodCfg.initialStateProjector,'function_handle')
                    yaPert = methodCfg.initialStateProjector(yaPert);
                end
                timeLeftGlobal = max(0, timeLimit - toc(startClock));
                attemptLimit = min(attemptLimitDefault, timeLeftGlobal);
                if attemptLimit <= 0
                    break
                end
                attemptStart = tic;
                try
                    [trajPert, infoPert] = integrateBdf2Trajectory(odeFun, domainMin, domainMax, yaPert, methodCfg, attemptStart, attemptLimit, utils);
                catch err
                    if strcmp(err.identifier,'displayBaseFn:timeout')
                        diagnostics.errorMessage = err.message;
                        diagnostics.warningId = err.identifier;
                        diagnostics.warningMsg = err.message;
                        sol = bestSol;
                        return
                    else
                        rethrow(err);
                    end
                end
                if infoPert.success && ~isempty(trajPert.x)
                    ybPert = trajPert.y(:,end);
                    rPertFull = bcFun(yaPert, ybPert);
                    totalBcEvals = totalBcEvals + 1;
                    if isempty(methodCfg.activeResidualIndices)
                        rPert = rPertFull;
                    else
                        rPert = rPertFull(methodCfg.activeResidualIndices);
                    end
                    if all(isfinite(rPert))
                        successCol = true;
                        break
                    end
                end
                delta = delta / 2;
            end
            if successCol
                J(:,col) = (rPert - residual) / delta;
            else
                J(:,col) = 0;
            end
        end
        if isempty(J)
            stepDir = zeros(numel(freeIdx),1);
        else
            if size(J,1) == size(J,2)
                rc = rcond(J);
                if ~isfinite(rc) || rc < 1e-12
                    stepDir = -pinv(J) * residual;
                else
                    stepDir = -J \ residual;
                end
            else
                stepDir = -pinv(J) * residual;
            end
        end
        ya(freeIdx) = ya(freeIdx) + stepDir;
        if ~isempty(methodCfg.initialStateProjector) && isa(methodCfg.initialStateProjector,'function_handle')
            ya = methodCfg.initialStateProjector(ya);
        end
    end

    if ~isempty(bestSol) && isfinite(bestResidual) && bestResidual <= shootTol
        sol = bestSol;
        diagnostics.status = 'success';
        diagnostics.errorMessage = sprintf('BDF2 shooting residual %.3e.', bestResidual);
        diagnostics.meshPoints = utils.getfieldWithDefault(bestInfo,'nodeCount',NaN);
        diagnostics.avgMeshStep = utils.getfieldWithDefault(bestInfo,'avgStep',NaN);
        diagnostics.iterations = utils.getfieldWithDefault(bestInfo,'totalNewtonIters',NaN);
        diagnostics.maxResidual = bestResidual;
        diagnostics.consoleLog = '';
        return
    end

    warning('displayBaseFn:bdf2Failure', ...
        'BDF2 solver failed to converge for [%g,%g]; skipping point (best residual %.3e vs tol %.3e).', ...
        domainMin, domainMax, bestResidual, shootTol);
    sol = [];
    diagnostics.status = 'fail';
    diagnostics.warningId = 'displayBaseFn:bdf2Failure';
    diagnostics.warningMsg = 'BDF2 solver failed to converge within limits.';
    diagnostics.errorMessage = diagnostics.warningMsg;
    diagnostics.maxResidual = bestResidual;
    diagnostics.consoleLog = '';
end

function vec = ensureColumnGuess(guessFcn, coord)
    val = guessFcn(coord);
    if isempty(val)
        error('displayBaseFn:guessFailure','Initial guess is empty at %g.', coord);
    end
    vec = val(:);
end

function [traj, info] = integrateBdf2Trajectory(odeFun, domainMin, domainMax, ya, methodCfg, startClock, timeLimit, utils)
    ya = ya(:);
    traj = struct('x', domainMin, 'y', ya, 'stats', struct('method','bdf2','nSteps',0,'nRejectedSteps',0,'nFevals',0));
    info = struct('success', false, 'errorMessage', '', 'nodeCount', 1, ...
                  'avgStep', NaN, 'totalNewtonIters', 0, 'acceptedSteps', 0, ...
                  'rejectedSteps', 0, 'fevals', 0, 'bcEvals', 0, ...
                  'relTolUsed', methodCfg.relativeTolerance, ...
                  'absTolUsed', methodCfg.absoluteTolerance, ...
                  'relaxations', 0);
    if domainMax <= domainMin
        info.success = true;
        return
    end

    relTol = methodCfg.relativeTolerance;
    absTol = methodCfg.absoluteTolerance;
    newtonTol = methodCfg.newtonTolerance;
    maxNewton = methodCfg.newtonMaxIterations;
    jacDelta = methodCfg.jacobianPerturbation;
    h = methodCfg.initialStepSize;
    hMin = methodCfg.minStepSize;
    hMax = methodCfg.maxStepSize;
    safety = max(0.1, min(0.95, methodCfg.stepSafetyFactor));
    growLim = max(1.1, methodCfg.stepGrowthLimit);
    shrinkLim = min(0.9, methodCfg.stepShrinkLimit);
    maxSteps = methodCfg.maxStepCount;

    odeEvalCount = 0;
    function fVal = evalOdeLocal(t, state)
        utils.enforceTimeLimit(@(~)0, startClock, timeLimit); %#ok<NASGU>
        odeEvalCount = odeEvalCount + 1;
        fVal = odeFun(t, state);
    end

    xCurr = domainMin;
    yCurr = ya;
    yPrev = [];
    fCurr = evalOdeLocal(xCurr, yCurr);

    xVals = xCurr;
    yVals = yCurr;
    accepted = 0;
    rejected = 0;
    totalNewton = 0;
    hasHistory = false;
    forceBootstrap = true;

    while xCurr < domainMax
        utils.enforceTimeLimit(@(~)0, startClock, timeLimit);
        if accepted + rejected >= maxSteps
            info.errorMessage = 'Maximum BDF2 step count exceeded.';
            break
        end

        h = min(max(h, hMin), min(hMax, domainMax - xCurr));
        if isfield(methodCfg,'frontierRefineLimit') && isfield(methodCfg,'frontierRefineStep') ...
                && ~isempty(methodCfg.frontierRefineLimit) && ~isempty(methodCfg.frontierRefineStep) ...
                && xCurr < methodCfg.frontierRefineLimit
            h = min(h, max(hMin, methodCfg.frontierRefineStep));
        end
        if h <= hMin * 0.5
            info.errorMessage = 'Step size underflow in BDF2 integrator.';
            break
        end

        if forceBootstrap || ~hasHistory
            [yCand, fCand, newtonIters, successStep] = solveBackwardEulerStep( ...
                xCurr, yCurr, fCurr, h, relTol, absTol, newtonTol, maxNewton, jacDelta, @evalOdeLocal);
            stepType = 'be';
        else
            [yCand, fCand, newtonIters, successStep] = solveBdf2Step( ...
                xCurr, yCurr, yPrev, h, fCurr, relTol, absTol, newtonTol, maxNewton, jacDelta, @evalOdeLocal);
            stepType = 'bdf2';
        end
        totalNewton = totalNewton + newtonIters;

        if ~successStep
            h = max(h * shrinkLim, hMin);
            forceBootstrap = true;
            rejected = rejected + 1;
            continue
        end

        switch stepType
            case 'be'
                errMetric = estimateBackwardEulerError(yCurr, yCand, fCand, h, absTol, relTol);
            otherwise
                errMetric = estimateBdf2Error(yCurr, yPrev, yCand, fCand, h, absTol, relTol);
        end
        if ~isfinite(errMetric)
            errMetric = Inf;
        end

        if errMetric > 1
            h = max(h * max(shrinkLim, safety * errMetric^(-0.5)), hMin);
            forceBootstrap = true;
            rejected = rejected + 1;
            continue
        end

        accepted = accepted + 1;
        xCurr = xCurr + h;
        xVals(end+1,1) = xCurr; %#ok<AGROW>
        yVals(:,end+1) = yCand; %#ok<AGROW>
        yPrev = yCurr;
        yCurr = yCand;
        fCurr = fCand;
        hasHistory = ~isempty(yPrev);
        forceBootstrap = false;

        if errMetric < 0.3
            proposed = min(h * min(growLim, max(1.1, safety * errMetric^(-0.5))), hMax);
        elseif errMetric > 0.8
            proposed = max(h * max(shrinkLim, safety * errMetric^(-0.5)), hMin);
        else
            proposed = h;
        end
        if abs(proposed - h) / max(h, eps) > 0.05
            forceBootstrap = true;
        end
        h = proposed;
    end

    traj.x = xVals.';
    traj.y = yVals;
    traj.stats = struct('method','bdf2', ...
                        'nSteps', accepted, ...
                        'nRejectedSteps', rejected, ...
                        'nFevals', odeEvalCount);
    info.success = abs(xCurr - domainMax) <= max(hMin, 1e-9);
    info.nodeCount = numel(xVals);
    info.avgStep = utils.averageMeshSpacing(xVals);
    info.acceptedSteps = accepted;
    info.rejectedSteps = rejected;
    info.totalNewtonIters = totalNewton;
    info.fevals = odeEvalCount;
    if ~info.success && isempty(info.errorMessage)
        info.errorMessage = 'BDF2 integrator halted before reaching the domain end.';
    end
end

function [yNext, fNext, newtonSteps, success] = solveBackwardEulerStep(x0, y0, f0, h, relTol, absTol, newtonTol, maxNewton, jacDelta, odeEvaluator)
    targetX = x0 + h;
    yGuess = y0 + h * f0;
    stateDim = numel(y0);
    success = false;
    fGuess = f0;

    for k = 1:maxNewton
        fGuess = odeEvaluator(targetX, yGuess);
        residual = yGuess - y0 - h * fGuess;
        if norm(residual, inf) <= newtonTol * (absTol + relTol * max(1, norm(yGuess, inf)))
            success = true;
            break
        end
        jac = approximateJacobian(targetX, yGuess, fGuess, jacDelta, odeEvaluator);
        delta = solveLinearSystem(eye(stateDim) - h * jac, residual);
        yGuess = yGuess - delta;
        if norm(delta, inf) <= newtonTol * max(1, norm(yGuess, inf))
            success = true;
            break
        end
    end
    newtonSteps = k;
    fNext = fGuess;
    yNext = yGuess;
end

function [yNext, fNext, newtonSteps, success] = solveBdf2Step(x0, y0, yPrev, h, f0, relTol, absTol, newtonTol, maxNewton, jacDelta, odeEvaluator)
    targetX = x0 + h;
    stateDim = numel(y0);
    yGuess = y0 + h * f0;
    success = false;
    fGuess = f0;

    for k = 1:maxNewton
        fGuess = odeEvaluator(targetX, yGuess);
        residual = (3*yGuess - 4*y0 + yPrev) / (2*h) - fGuess;
        if norm(residual, inf) <= newtonTol * (absTol + relTol * max(1, norm(yGuess, inf)))
            success = true;
            break
        end
        jac = approximateJacobian(targetX, yGuess, fGuess, jacDelta, odeEvaluator);
        delta = solveLinearSystem((3/(2*h))*eye(stateDim) - jac, residual);
        yGuess = yGuess - delta;
        if norm(delta, inf) <= newtonTol * max(1, norm(yGuess, inf))
            success = true;
            break
        end
    end

    newtonSteps = k;
    fNext = fGuess;
    yNext = yGuess;
end

function errMetric = estimateBackwardEulerError(yCurr, yNext, fNext, h, absTol, relTol)
    defect = yNext - yCurr - h * fNext;
    scale = absTol + relTol * max(norm(yNext, inf), norm(yCurr, inf));
    if scale == 0
        scale = 1;
    end
    errMetric = norm(defect, inf) / scale;
end

function errMetric = estimateBdf2Error(yCurr, yPrev, yNext, fNext, h, absTol, relTol)
    predictor = yCurr + h * fNext;
    errVec = yNext - predictor;
    scale = absTol + relTol * max([norm(yNext, inf), norm(yCurr, inf), norm(yPrev, inf)]);
    if scale == 0
        scale = 1;
    end
    errMetric = norm(errVec, inf) / scale;
end

function jac = approximateJacobian(x, y, fBase, delta, odeEvaluator)
    n = numel(y);
    jac = zeros(n);
    for idx = 1:n
        pert = delta * max(1, abs(y(idx)));
        if pert == 0
            pert = delta;
        end
        yPert = y;
        yPert(idx) = yPert(idx) + pert;
        fPert = odeEvaluator(x, yPert);
        jac(:, idx) = (fPert - fBase) / pert;
    end
end

function delta = solveLinearSystem(A, rhs)
    if any(isnan(rhs))
        delta = rhs;
        return
    end
    rc = rcond(A);
    if ~isfinite(rc) || rc < 1e-12
        delta = pinv(A) * rhs;
    else
        delta = A \ rhs;
    end
end

function ya = seedInitialStateWithCollocation(ya, odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel, utils)
    rel = utils.getfieldWithDefault(methodCfg,'seedCollocationRelTol',methodCfg.relativeTolerance);
    absTol = utils.getfieldWithDefault(methodCfg,'seedCollocationAbsTol',methodCfg.absoluteTolerance);
    maxPts = utils.getfieldWithDefault(methodCfg,'seedCollocationMaxPoints',methodCfg.maxMeshPoints);
    seedCfg = methodCfg;
    seedCfg.relativeTolerance = rel;
    seedCfg.absoluteTolerance = absTol;
    seedCfg.maxMeshPoints = maxPts;
    seedCfg.maxSolverAttempts = 2;
    seedCfg.timeLimitSeconds = min(5, methodCfg.timeLimitSeconds/2);
    try
        [seedSol, ~] = collocationSolverBase(odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, seedCfg, domainLabel, 'bvp4c', utils);
        if ~isempty(seedSol) && isfield(seedSol,'y') && ~isempty(seedSol.y)
            ya = seedSol.y(:,1);
        end
    catch err
        warning('displayBaseFn:bdf2Seed','Collocation seed failed: %s', err.message);
    end
end
