function [sol, diagnostics] = bvp5cSolver(odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel, utils)
% Objective: Wrapper exposing bvp5c through the shared solver interface.
% Purpose : Keep the dispatcher generic while reusing collocation infrastructure.
% SWOT     S: minimal code; W: inherits bvp5c limitations; O: future tweaks centralized; T: solver removal affects callers.

    [sol, diagnostics] = collocationSolverBase( ...
        odeFun, bcFun, guessFcn, domainMin, domainMax, stepSize, methodCfg, domainLabel, 'bvp5c', utils);
end
