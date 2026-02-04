function values = values()
%VALUES Parameter definition for Arif's configuration.
% figure 3 might be terbalik 260102_1517_arif_bvp5c_m_0p1v0p3 but figure 4
% is right
% sweep: (1) change inclination angle pi/3 (tak dpt 2nd sol), pi/4, pi/6 (Anuar)


% Hybrid Nanofluid Al-Cu with based fluid EG

% Hybrid Nanofluid Al-Cu with based fluid EG
p = struct();
p.Bi = 0.2; % Biot Number
p.A = -1.2; % Unsteady parameter
p.M = 0.1; % Magnetic Field Parameter
p.lambda = 0.5; % bouyancy parameter
p.Rd = 0.2; % Radiation Parameter
p.S = 2; % Suction/Injection Parameter
p.Sl = 0.3; % Velocity Slip parameter
p.alpha = pi/4; % Inclination Angle

p.delta = 0.2; % Material parameter, Thickness of Thermal boundary layer
p.omega = 0.1; % Material parameter, Dimensionless heat transfer factor
p.Pr = 204; % Prandtl Number, EG % Aminuddin impact of ...

m = struct();
m.phi1 = 0.01; % Volume fraction for Alumina
m.phi2 = 0.01; % Volume fraction for Copper
m.phiHnf = m.phi1 + m.phi2; % Hybrid volume fraction
m.rhoS1 = 3970; % Density of Alumina
m.rhoS2 = 8933; % Density of Copper
m.rhoF = 1114; % Density of EG
m.betaS1 = 0.85e-5; % Thermal Expension of Al
m.betaS2 = 1.67e-5; % Thermal Expension of Cu
m.betaF = 57e-5; % Density of EG
% m.sigmaS1 = 3.5e7; % Electrical conductivity of Al metal
m.sigmaS1 = 0; % Electrical conductivity of Al
m.sigmaS2 = 5.96e7;  % Electrical conductivity of Cu
m.sigmaF = 5.5e-6;  % Electrical conductivity of EG
m.CpS1 = 765; % Heat capacity at constant pressure Al
m.CpS2 = 385; % Heat capacity at constant pressure Cu
m.CpF = 2415; % Heat capacity at constant pressure EG
m.kS1 = 40; % Thermal Conductivity Al
m.kS2 = 401; % Thermal Conductivity Cu
m.kF = 0.252; % Thermal Conductivity EG

n = deriveNFromM(m);

values = struct();
values.p = p;
values.m = m;
values.n = n;
values.deriveNFromM = @deriveNFromM;
values.guesses = {
    @(x) firstBranchGuess(x, p);
    @(x) secondBranchGuess(x, p);
    };
values.domainMinList = [0, 0];
values.domainMaxList = [10, 10];

values.sweep = struct();
values.sweep.mValues = [0.2];
% lambdaVals=linspace(-1.85,1, 100);
lambdaVals=linspace(-2,1, 4);


lambdaVals1 = linspace(0.19827586,1, 5);

lambdaVals2 = linspace(0.04958789,0.19827586, 5);
lambdaVals3 = linspace(0.03267573,0.04958789, 5);
lambdaVals4 = linspace(0.02703834,0.03267573, 5);
lambdaVals5 = linspace(0.01211382,0.02703834, 5);

% lambdaVals5DeadLine = linspace(-0.00833887,0.01211382, 15);
lambdaVals6 = linspace(-0.00905697,-0.00833887, 5);
lambdaVals7 = linspace(-0.01121125,-0.00905697, 5);
lambdaVals8 = linspace(-0.01395442,-0.01121125, 5);
lambdaVals9 = linspace(-0.01669759,-0.01395442, 5);
lambdaVals10 = linspace(-0.14916112,-0.01669759, 5);

lambdaVals11 = linspace(-0.28162466,-0.14916112, 5);

lambdaVals12 = linspace(-0.54655172,-0.28162466, 5);
lambdaVals13 = linspace(-0.73275862,-0.54655172, 5);
lambdaVals14 = linspace(-0.86093735,-0.73275862, 5);

lambdaVals15 = linspace(-1.01206897,-0.86093735, 5);
lambdaVals16 = linspace(-1.19827586,-1.01206897, 5);
lambdaVals17 = linspace(-1.30723845,-1.19827586, 5);

lambdaVals18 = linspace(-1.57068966,-1.30723845, 5);
lambdaVals19 = linspace(-1.75689655,-1.57068966, 5);
lambdaVals20 = linspace(-1.82671604,-1.75689655, 5);

lambdaVals21 = linspace(-1.85,-1.82671604, 5);

% the critical point atleast converge to the roundoff 7 decimal places
% need to verify critical point figure 4
% lambdaVals = unique([
%     % lambdaVals(:)
%     % lambdaVals1(:);
%     % lambdaVals2(:);
%     % lambdaVals3(:);
%     % lambdaVals4(:);
%     % % lambdaVals5DeadLine(:);
%     % lambdaVals6(:);
%     % lambdaVals7(:);
%     % lambdaVals8(:);
%     % lambdaVals9(:);
%     % lambdaVals10(:);
%     % lambdaVals11(:);
%     % lambdaVals12(:);
%     % lambdaVals13(:);
%     % lambdaVals14(:);
%     % lambdaVals15(:);
%     % lambdaVals16(:);
%     % lambdaVals17(:);
%     % lambdaVals18(:);
%     % lambdaVals19(:);
%     % lambdaVals20(:);
%     % lambdaVals21(:);
%     % % -1.82671604;
%     % % -1.30723845;
%     % % -0.86093735;
%     % % -0.01669759;
%     % % -0.01121125;
%     % % -0.00833887;
%     % % 0.01211382;
%     % % 0.02703834;
%     % % 0.04958789;
%     % 0
%     ], 'sorted');
% (lambdaVals >= min & lambdaVals <= max)
% values.sweep.excludeLambdaRanges = [
%     -0.01669759 -0.00833887
%     -1.6 -1.4
%     ];
values.sweep.excludeLambdaRanges = [];
values.sweep.lambdaVals = lambdaVals;
values.sweep.primary = struct( ...
    'paramScope', 'p', ...
    'paramName', 'M', ...
    'name', 'M', ...
    'values', values.sweep.mValues, ...
    'labelFcn', @(val) sprintf('M=%.2f', val));
values.sweep.secondary = struct( ...
    'paramScope', 'p', ...
    'paramName', 'A', ...
    'name', 'Unsteady, A', ...
    'values', values.sweep.lambdaVals, ...
    'labelFcn', @(val) sprintf('\\A=%.2f', val), ...
    'probeSides', true, ...
    'maxConsecutiveFails', 1, ...
    'refine', struct('onFail', true, 'numSubdiv', 5), ...
    'useSmartLinspace', true, ...
    'excludeRanges', values.sweep.excludeLambdaRanges, ...
    'group', struct( ...
    'paramScope', 'p', ...
    'paramName', 'M', ...
    'values', values.sweep.mValues, ...
    'labelFcn', @(val) sprintf('M=%.2f', val)));
% values.sweep.smartLinspaceLambdaVals = smartLinspace(-1.85, 1, 100, 0.08, 0.1);
end

function n = deriveNFromM(m)
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
end

function guess = firstBranchGuess(x, p)
x = x(:).';
decay = exp(-x);
lambdaShift = 0.04 * tanh(-p.lambda);
profileShift = exp(-0.5 * x);
f = p.S + (1 - decay);
fp = decay;
fpp = -decay;
shift = lambdaShift * profileShift .* (1 + 0.5 * x);
f = f + shift;
fp = fp + lambdaShift * profileShift .* (0.5 - 0.25 * x);
fpp = fpp + lambdaShift * profileShift .* (-0.5 + 0.125 * x);
thetaBase = decay;
thetaCorrection = lambdaShift * profileShift;
guess = [
    f;
    fp;
    fpp;
    thetaBase + thetaCorrection;
    -thetaBase + thetaCorrection
    ];
end

function guess = secondBranchGuess(x, p)
x = x(:).';
fast = exp(-x);
slow = exp(-0.1 * x);
fp = 1 + 3 * slow - 6 * fast;
fpp = -0.3 * slow + 6 * fast;
f = p.S + x + 30 * (1 - slow) - 6 * (1 - fast);
thetaBase = -exp(-0.1 * x);
thetaCorrection = 0.1 * exp(-0.1 * x);
guess = [
    f;
    fp;
    fpp;
    thetaBase;
    thetaCorrection
    ];
end
