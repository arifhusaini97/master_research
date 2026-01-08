function spec = method(cfg)
%METHOD Numerical strategy for Arif profile (bvp4c collocation).

if nargin < 1
    cfg = struct();
end

spec = struct();
spec.solver = 'bvp5c';
spec.displayName = 'Collocation (bvp5c)';
spec.relativeTolerance = 1e-10;
spec.absoluteTolerance = 1e-10;
% default like dr_zana template is minMeshPoints = 5
spec.minMeshPoints = 51;
spec.domainGridSizeList = [51, 51];
spec.maxMeshPoints = 60000;
end
