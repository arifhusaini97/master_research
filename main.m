function summary = main(problemName)
%MAIN Entry point: select profile, run sweeps, emit artifacts via one call.

if nargin < 1 || strlength(string(problemName)) == 0
  problemName = 'arif';
end

format long g

cfg = factory('loadProblemValues', problemName);
model = factory('loadProblemModel', problemName);
method = factory('loadProblemMethod', problemName, cfg);

summary = factory('runProblem', cfg, model, method);

fprintf('Problem "%s" complete. Artifacts stored in %s\n', ...
  summary.problem, summary.outputDir);
end
