function [sol, diagnostics] = bvp4cSolver(odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel, utils)
% Objective: Wrapper exposing bvp4c through the shared solver interface.
% Purpose : Keep displayBaseFn agnostic from specific MATLAB solvers.
% SWOT     S: plug-and-play; W: still dependent on MATLAB's bvp4c; O: easy to tweak options; T: any bvp4c bug propagates.

    [sol, diagnostics] = collocationSolverBase( ...
        odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel, 'bvp4c', utils);
end
