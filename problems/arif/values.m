function cfg = values()
%VALUES Parameter definition for Arif's configuration.

cfg = struct();
p = struct();
p.delta = 0.2;
p.omega = 0.1;
p.A = -1.2;
p.M = 0.01;
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

branchOneState = [0.1; 0.816353; -0.612155; -0.126812; -0.225362];
branchTwoState = [0.1; 0.93423; -0.219232; 1.53037; 0.106074];

branchSpecs = [
    struct('state', branchOneState, 'fpDecay', 0.9, 'fpTail', 0.01, ...
    'fppDecay', 0.7, 'fppTail', -0.02, 'tFar', -0.05, 'tDecay', 0.55);
    struct('state', branchTwoState, 'fpDecay', 0.6, 'fpTail', 0.02, ...
    'fppDecay', 0.45, 'fppTail', -0.05, 'tFar', 0.2, 'tDecay', 0.4)
    ];

cfg.guesses = {
    @(x) branchGuessProfile(x, branchSpecs(1), p);
    @(x) branchGuessProfile(x, branchSpecs(2), p);
    };

cfg.domainMinList = [0, 0];
cfg.domainMaxList = [3.5, 3.5];
cfg.domainGridSizeList = [5, 5];

cfg.p = p;
cfg.m = m;
cfg.n = n;
end

function profile = branchGuessProfile(x, spec, p)
eta = x(:).';
if isempty(eta)
    eta = 0;
end
state = spec.state(:);
if numel(state) ~= 5
    error('values:invalidSeed','State vector must have 5 components.');
end
fpDecay = max(spec.fpDecay, 0.1);
fpTail = spec.fpTail;
fppDecay = max(spec.fppDecay, 0.1);
fppTail = spec.fppTail;
tFar = spec.tFar;
tDecay = max(spec.tDecay, 0.2);

fpp = fppTail + (state(3) - fppTail) .* exp(-fppDecay .* eta);
fpCore = fpTail + (state(2) - fpTail) .* exp(-fpDecay .* eta);
fp = 1 + p.Sl * fpp;
fp = 0.7 * fp + 0.3 * fpCore;
fp(1) = state(2);

if fpDecay <= 0
    f = p.S + fpTail .* eta;
else
    f = p.S + fpTail .* eta + ((state(2) - fpTail) / fpDecay) .* (1 - exp(-fpDecay .* eta));
end

t = tFar + (state(4) - tFar) .* exp(-tDecay .* eta);
tp = -tDecay * (state(4) - tFar) .* exp(-tDecay .* eta);

profile = [
    f;
    fp;
    fpp;
    t;
    tp
    ];
end
