function values = values()
%VALUES Parameter definition for Arif's configuration.
% figure 3 might be terbalik 260102_1517_arif_bvp5c_m_0p1v0p3 but figure 4
% is right
% sweep: (1) change inclination angle pi/3 (tak dpt 2nd sol), pi/4, pi/6 (Anuar)

values = struct();
p = struct();
p.delta = 0.2;
p.omega = 0.1;
p.A = -1.2;
p.M = 0.1;
p.lambda = -0.2;
p.Pr = 6.2;
p.Rd = 0.2;
p.S = 0.1;
p.Sl = 0.3;
p.Bi = 0.2;
p.alpha = pi/4;

m = struct();
m.phi1 = 0.01; m.phi2 = 0.01;
m.phiHnf = m.phi1 + m.phi2;
m.rhoS1 = 3970; m.rhoS2 = 8933; m.rhoF = 997.1;
m.betaS1 = 0.85e-5; m.betaS2 = 1.67e-5; m.betaF = 21;
m.sigmaS1 = 3.5e7; m.sigmaS2 = 5.96e7; m.sigmaF = 0.05;
m.CpS1 = 765; m.CpS2 = 385; m.CpF = 4179;
m.kS1 = 40; m.kS2 = 401; m.kF = 0.613;

n = struct();
n.bMu = 1 / ((1 - m.phiHnf)^2.5);
n.bRho = (1 - m.phiHnf) + m.phi1*(m.rhoS1/m.rhoF) + m.phi2*(m.rhoS2/m.rhoF);
n.bBeta = ((1 - m.phiHnf)*m.rhoF + m.phi1*m.rhoS1*(m.betaS1/m.betaF) + m.phi2*m.rhoS2*(m.betaS2/m.betaF)) ...
    * (1 / (n.bRho * m.rhoF));
n.bSigma = (((m.phi1*m.sigmaS1 + m.phi2*m.sigmaS2)/m.phiHnf) + 2*m.sigmaF + 2*(m.phi1*m.sigmaS1 + m.phi2*m.sigmaS2) - 2*m.phiHnf*m.sigmaF) ...
    / (((m.phi1*m.sigmaS1 + m.phi2*m.sigmaS2)/m.phiHnf) + 2*m.sigmaF - (m.phi1*m.sigmaS1 + m.phi2*m.sigmaS2) + m.phiHnf*m.sigmaF);
n.bRhoCp = (1 - m.phiHnf) + m.phi1*(m.rhoS1*m.CpS1)/(m.rhoF*m.CpF) + m.phi2*(m.rhoS2*m.CpS2)/(m.rhoF*m.CpF);
n.bK = (((m.phi1*m.kS1 + m.phi2*m.kS2)/m.phiHnf) + 2*m.kF + 2*(m.phi1*m.kS1 + m.phi2*m.kS2) - 2*m.phiHnf*m.kF) ...
    / (((m.phi1*m.kS1 + m.phi2*m.kS2)/m.phiHnf) + 2*m.kF - (m.phi1*m.kS1 + m.phi2*m.kS2) + m.phiHnf*m.kF);

values.p = p;
values.m = m;
values.n = n;
values.guesses = {
    @(x) firstBranchGuess(x, p);
    @(x) secondBranchGuess(x, p);
    };
values.domainMinList = [0, 0];
values.domainMaxList = [10, 10];

values.sweep = struct();
values.sweep.mValues = [0.1,0.2,0.3];
lambdaVals = linspace(-1.85, 0.85, 100);
% the critical point atleast converge to the roundoff 7 decimal places
lambdaVals = unique([lambdaVals, -1.8268585, -1.30733462342,-0.860944740656, -0.01669759,-0.01121125,-0.00833887,0.01211382,0.02703834, 0.04958789], 'sorted');
% (lambdaVals >= min & lambdaVals <= max)
% values.sweep.excludeLambdaRanges = [
%     -0.01669759 -0.00833887
%     -1.6 -1.4
%     ];
values.sweep.excludeLambdaRanges = [];
values.sweep.lambdaVals = lambdaVals;
end

function guess = firstBranchGuess(x, p)
x = x(:).';
decay = exp(-x);
guess = [
    p.S + 0.2 + decay;
    decay;
    -decay;
    decay;
    decay
    ];
end

function guess = secondBranchGuess(x, p)
x = x(:).';
decay = exp(-x);
guess = [
    p.S + (1 - decay);
    decay;
    -decay;
    decay;
    -decay
    ];
end
