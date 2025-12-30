function [sol, diagnostics] = collocationSolverBase(odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel, solverName, utils)
% Objective: Shared driver for MATLAB's collocation solvers (bvp4c/bvp5c).
% Purpose : Normalize retry/time-limit logic once so individual wrappers stay thin.
% SWOT     S: single implementation; W: relies on MATLAB solvers availability;
%          O: easy to add new collocation variants; T: regressions affect all callers.

    diagnostics = utils.initBranchDiagnostics(stepSize);
    mesh = linspace(domainMin, domainMax, max(stepSize, 51));
    rel = methodCfg.relativeTolerance;
    absTol = methodCfg.absoluteTolerance;
    maxMeshPts = methodCfg.maxMeshPoints;
    opts = bvpset('Stats','on','RelTol',rel,'AbsTol',absTol,'NMax',maxMeshPts);
    maxAttempts = methodCfg.maxSolverAttempts;
    timeoutSeconds = methodCfg.timeLimitSeconds;

    startClock = tic;
    timedOde = @(x,y) utils.enforceTimeLimit(odeFun, startClock, timeoutSeconds, x, y);
    timedBc  = @(ya,yb) utils.enforceTimeLimit(bcFun, startClock, timeoutSeconds, ya, yb);

    warnState = warning;
    warning('off','all');
    warnCleanup = onCleanup(@() warning(warnState)); %#ok<NASGU>

    solverLabel = upper(solverName);
    lastError = [];

    for attempt = 1:maxAttempts
        diagnostics.attempts = attempt;
        if toc(startClock) > timeoutSeconds
            lastError = MException('displayBaseFn:timeout', ...
                                   'Time limit of %.1f s exceeded before attempt %d.', timeoutSeconds, attempt);
            break;
        end

        guessMatrix = [];
        try
            guessMatrix = utils.sampleInitialGuess(guessFcn, mesh);
        catch
            guessMatrix = [];
        end
        if isempty(mesh)
            diagnostics.initialGuessMesh = [];
        else
            diagnostics.initialGuessMesh = mesh(1);
        end
        if ~isempty(guessMatrix)
            diagnostics.initialGuessProfile = guessMatrix(:,1);
        else
            diagnostics.initialGuessProfile = [];
        end
        solinit = bvpinit(mesh, guessFcn);
        try
            lastwarn('', '');
            [~, solCandidate] = evalc(sprintf('%s(timedOde, timedBc, solinit, opts)', solverName));
            [warnMsg, warnId] = lastwarn;
            lastwarn('', '');

            if isCollocationWarning(warnId, solverName)
                lastError = MException(warnId, warnMsg);
                mesh = sort([mesh, utils.midpoints(mesh)]);
                rel = min(1e-5, rel*5);
                absTol = min(1e-7, absTol*5);
                opts = bvpset(opts,'RelTol',rel,'AbsTol',absTol,'NMax',maxMeshPts);
                diagnostics.warningId = warnId;
                diagnostics.warningMsg = warnMsg;
                continue;
            end

            stats = struct();
            if isfield(solCandidate,'stats')
                stats = solCandidate.stats;
            end
            meshPoints = numel(solCandidate.x);
            diagnostics.status = 'success';
            diagnostics.warningId = warnId;
            diagnostics.warningMsg = warnMsg;
            diagnostics.meshPoints = meshPoints;
            diagnostics.avgMeshStep = utils.averageMeshSpacing(solCandidate.x);
            diagnostics.maxResidual = utils.extractMaxResidual(solCandidate);
            diagnostics.errorMessage = '';
            diagnostics.initialGuessError = utils.computeInitialGuessError(solCandidate, guessFcn);
            diagnostics.iterations = utils.extractIterationCount(solCandidate);

            odeEvalCount = extractStatCount(stats, {'nODEevals','nOdeEvals','nfevals','nFevals'});
            bcEvalCount = extractStatCount(stats, {'nBCevals','nBcEvals'});
            logBody = utils.buildMethodConsoleLog(meshPoints, diagnostics.maxResidual, odeEvalCount, bcEvalCount);
            meta = methodCfg;
            meta.solutionRange = [solCandidate.x(1), solCandidate.x(end)];
            diagnostics.consoleLog = utils.formatSolverConsoleLog(logBody, attempt, domainMin, domainMax, domainLabel, meta);
            if ~isempty(diagnostics.consoleLog)
                utils.emitRunLog('%s\n', diagnostics.consoleLog);
            end
            sol = solCandidate;
            return;
        catch err
            if strcmp(err.identifier,'displayBaseFn:timeout')
                lastError = err;
                break;
            end
            lastError = err;
            mesh = sort([mesh, utils.midpoints(mesh)]);
            rel = min(1e-5, rel*5);
            absTol = min(1e-7, absTol*5);
            opts = bvpset(opts,'RelTol',rel,'AbsTol',absTol,'NMax',maxMeshPts);
            diagnostics.warningId = err.identifier;
            diagnostics.warningMsg = err.message;
            lastwarn('', '');
        end
    end

    if isempty(lastError)
        lastError = MException('displayBaseFn:solveBranchFail','Unknown failure.');
    end
    warning('displayBaseFn:solveBranchFail', ...
        '%s failed after %d attempts: %s', solverLabel, maxAttempts, lastError.message);
    sol = [];
    diagnostics.status = 'fail';
    diagnostics.errorMessage = lastError.message;
end

function count = extractStatCount(stats, candidates)
    count = NaN;
    if isempty(stats) || ~isstruct(stats)
        return
    end
    for k = 1:numel(candidates)
        fld = candidates{k};
        if isfield(stats, fld) && ~isempty(stats.(fld))
            count = stats.(fld);
            return
        end
    end
end

function flag = isCollocationWarning(warnId, solverName)
% Objective: Identify tolerable warnings emitted by MATLAB's collocation solvers.
% Purpose : Allows adaptive retries when solver reports mesh issues instead of aborting immediately.
% SWOT     S: isolates warning logic; W: relies on message text; O: tweak for new solvers; T: new warn IDs may slip.
    if nargin < 2 || strlength(string(solverName)) == 0
        solverName = 'bvp4c';
    end
    flag = ~isempty(warnId) && contains(lower(warnId), lower(string(solverName)));
end
