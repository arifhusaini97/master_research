function fns = model()
%MODEL Encapsulates Natasha-specific solver functions.

fns = struct();
fns.ode = @natasha_odeFun;
fns.bc  = @natasha_bcFun;
fns.localSkinFriction = @natasha_localSkinFriction;
fns.nusseltNumber = @natasha_nusseltNumberFn;
end

function ff = natasha_odeFun(x, y, p, n)
f   = y(1);
fp  = y(2);
fpp = y(3);
th  = y(4);
thp = y(5);

denom = n.muRatio + p.B - p.E * p.B * (fpp^2);
numer = n.sigmaRatio * p.M * fp + n.rhoRatio * (fp^2 - f * fpp);

ff = [
    fp
    fpp
    numer / denom
    thp
    natasha_thermalRHS(f,fp,th,thp,x,p,n)
    ];
end

function q = natasha_thermalRHS(f,fp,t,tp,x,p,n)
unsteady = 0;
if isfield(p,'A') && ~isempty(p.A)
    unsteady = p.A * (2*t + 0.5*x*tp);
end
convective = fp * t - f * tp + unsteady;
q = p.Pr * n.rhoCpRatio / n.kRatio * convective;
end

function res = natasha_bcFun(ya, yb, p, ~)
res = [
    ya(1)-p.S
    ya(2)-p.lambda-p.L1*((1-p.L2*ya(3))^(-1/2))*ya(3)
    ya(4)-1
    yb(2)
    yb(4)
    ];
end

function CfRe12 = natasha_localSkinFriction(sol, p, n)
y0    = factory('evalSolution', sol, 0);      % [f, f', f'', t, t']
fpp0  = y0(3);

muStar = n.muRatio + p.B;
powellCorr = p.B * p.E / 3;

CfRe12 = muStar * fpp0 - powellCorr * (fpp0^3);
end

function NuRe12 = natasha_nusseltNumberFn(sol, ~, n)
y0  = factory('evalSolution', sol, 0);      % [f, f', f'', t, t']
tp0 = y0(5);
NuRe12 = -n.kRatio * tp0;
end
