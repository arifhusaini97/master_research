function fns = model()
%MODEL Encapsulates Arif-specific solver functions.

fns = struct();
fns.ode = @arif_odeFun;
fns.bc  = @arif_bcFun;
fns.localSkinFriction = @arif_localSkinFriction;
fns.nusseltNumber = @arif_nusseltNumberFn;
end

function ff = arif_odeFun(x, y, p, n)
denom = (p.omega * (p.delta * y(3)^2 - 1)) - n.bMu;
inertial = y(2)^2 - y(1)*y(3) + p.A*(y(2) + (x/2)*y(3));
electromagnetic = n.bSigma * p.M * y(2);
rhsMomentum = -n.bRho*inertial + n.bRho*n.bBeta * p.lambda * y(4) * cos(p.alpha);

thermalNumer = y(2)*y(4) - y(1)*y(5) + p.A*(2*y(4) + (x/2)*y(5));
thermalDenom = n.bK + (4/3)*p.Rd;

ff = [
    y(2)
    y(3)
    (rhsMomentum - electromagnetic) / denom
    y(5)
    p.Pr * n.bRhoCp * thermalNumer / thermalDenom
    ];
end

function res = arif_bcFun(ya, yb, p, n)
% res = [
%     ya(1)-p.S
%     ya(2)-1-(p.Sl*ya(3))
%     n.bK-ya(5)+(p.Bi*(1-ya(4)))
%     yb(2)
%     yb(4)
%     ];


res = [ ya(1)-p.S
    ya(2)-1-p.Sl*ya(3)
    %ya(4)-1
    ya(5)+p.Bi*(1-ya(4))
    yb(2)
    yb(4)];
end

function CfRe12 = arif_localSkinFriction(sol, p, n)
y0    = factory('evalSolution', sol, 0);
fpp0  = y0(3);
CfRe12 = (n.bMu + p.omega)*fpp0 - (1/3)*p.delta*p.omega*(fpp0^3);
end

function NuRe12 = arif_nusseltNumberFn(sol, p, n)
y0  = factory('evalSolution', sol, 0);
tp0 = y0(5);
NuRe12 = -(n.bK + (4/3)*p.Rd)*tp0;
end
