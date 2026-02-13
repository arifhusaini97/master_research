function CodeA

format long g

%Define all parameters
global p m n

%Boundary layer thickness & stepsize
etaMin = 0;
etaMax1 = 10;
etaMax2 = 10;
stepsize1 = 51;
stepsize2 = 51;


% Hybrid Nanofluid Al-Cu with based fluid EG
p = struct();
p.Bi = 0.1; % Biot Number %Okello 0.1
p.A = 0.2; % Unsteady parameter %Aziz 0.2
p.M = 0.1; % Magnetic Field Parameter
p.lambda = -1; % bouyancy parameter %Anuar -1
p.Rd = 0.2; % Radiation Parameter % Aziz 0.2
p.S = 1; % Suction/Injection Parameter
p.Sl = 0.5; % Velocity Slip parameter
p.alpha = pi/4; % Inclination Angle

p.delta = 0.2; % Material parameter, Thickness of Thermal boundary layer
p.omega = 0.4; % Material parameter, Dimensionless heat transfer factor %Hahim 0.4
p.Pr = 204; % Prandtl Number, EG %drZana. 6.2 is for water

m = struct();
m.phi1 = 0.01; % Volume fraction for Alumina
m.phi2 = 0.01; % Volume fraction for Copper
m.phiHnf = m.phi1 + m.phi2; % Hybrid volume fraction
m.rhoS1 = 3970; % Density of Alumina
m.rhoS2 = 8933; % Density of Copper
m.rhoF = 1114; % Density of EG
m.betaS1 = 0.85e-5; % Thermal Expension of Al
m.betaS2 = 1.67e-5; % Thermal Expension of Cu
m.betaF = 5.7e-4; % Density of EG
m.sigmaS1 = 3.5e7; % Electrical conductivity of Al2O3
% m.sigmaS1 = 1e-10; % Electrical conductivity of Al2O3 insulator
m.sigmaS2 = 5.96e7;  % Electrical conductivity of Cu
m.sigmaF = 5.5e-6;  % Electrical conductivity of EG
m.CpS1 = 765; % Heat capacity at constant pressure Al
m.CpS2 = 385; % Heat capacity at constant pressure Cu
m.CpF = 2415; % Heat capacity at constant pressure EG
m.kS1 = 40; % Thermal Conductivity Al
m.kS2 = 401; % Thermal Conductivity Cu
m.kF = 0.252; % Thermal Conductivity EG

n = struct();

n.phiMu = 1 / ((1 - m.phiHnf)^2.5);
n.phiRho = (1 - m.phiHnf) + ((m.phi1*m.rhoS1*m.CpS1 + m.phi2*m.rhoS2*m.CpS2)/(m.rhoF *m.CpF));
n.rhoHnf=n.phiRho*m.rhoF;

n.phiBeta = (1/n.rhoHnf)*(((m.phi1*m.rhoS1*m.betaS1 + m.phi2*m.rhoS2*m.betaS2)/m.betaF)+((1-m.phiHnf)*m.rhoF));

n.phiSigma = (((m.phi1*m.sigmaS1 + m.phi2*m.sigmaS2)/m.phiHnf) + 2*m.sigmaF + 2*(m.phi1*m.sigmaS1 + m.phi2*m.sigmaS2) - 2*m.phiHnf*m.sigmaF) ...
    / (((m.phi1*m.sigmaS1 + m.phi2*m.sigmaS2)/m.phiHnf) + 2*m.sigmaF - (m.phi1*m.sigmaS1 + m.phi2*m.sigmaS2) + m.phiHnf*m.sigmaF);
n.phiRhoCp = (1-m.phiHnf)+(((m.phi1*m.rhoS1*m.CpS1)+(m.phi2*m.rhoS2*m.CpS2))/(m.rhoF*m.CpF))
n.phiK = (((m.phi1*m.kS1 + m.phi2*m.kS2)/m.phiHnf) + 2*m.kF + 2*(m.phi1*m.kS1 + m.phi2*m.kS2) - 2*m.phiHnf*m.kF) ...
    / (((m.phi1*m.kS1 + m.phi2*m.kS2)/m.phiHnf) + 2*m.kF - (m.phi1*m.kS1 + m.phi2*m.kS2) + m.phiHnf*m.kF);

%%%%%%%%%%%%%%%%%%%%%%   first solution   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options =  bvpset('stats','off','RelTol',1e-6,'AbsTol',1e-8);
solinit = bvpinit (linspace (etaMin, etaMax1, stepsize1), @OdeInit1);
sol = bvp5c (@OdeBVP, @OdeBC, solinit, options);
eta = linspace (etaMin, etaMax1, stepsize1);
y = deval (sol, eta);

figure(1)
plot(sol.x,sol.y(2,:),'r')
xlabel('\eta')
ylabel('f`(\eta)')
hold on

figure(2)
plot(sol.x,sol.y(4,:),'r')
xlabel('\eta')
ylabel('t(\eta)')
hold on

%Saving the output in txt file for first solution
descris = [sol.x; sol.y];
save 'upper.txt' descris -ascii

%Displaying the output for first solution
fprintf('\nFirst solution:\n');
fprintf('f"(0) = %7.9f\n',y(3));
fprintf('-t`(0) = %7.9f\n',-y(5));
CfRe12 = localSkinFriction(sol, p, n);
NuRe12 = nusseltNumberFn(sol, p, n);
fprintf('Cf Re^{1/2} = %7.9f\n',CfRe12);
fprintf('Nu Re^{-1/2} = %7.9f\n',NuRe12);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%   second solution   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = bvpset('stats','off','RelTol',1e-6,'AbsTol',1e-8,'NMax',5000);
solinit = bvpinit (linspace (etaMin, etaMax2, stepsize2), @OdeInit2);
sol2 = bvp5c (@OdeBVP, @OdeBC, solinit, options);
eta = linspace (etaMin, etaMax2, stepsize2);
y = deval (sol2, eta);

figure(1)
plot(sol2.x,sol2.y(2,:),'--r')
xlabel('\eta')
ylabel('f`(\eta)')
hold on

figure(2)
plot(sol2.x,sol2.y(4,:),'--r')
xlabel('\eta')
ylabel('t(\eta)')
hold on

%Saving the output in txt file for second solution
descris = [sol2.x; sol2.y];
save 'lower.txt' descris -ascii

%Displaying the output for second solution
fprintf('\nSecond solution:\n');
fprintf('f"(0) = %7.9f\n',y(3));
fprintf('-t`(0) = %7.9f\n',-y(5));
CfRe12 = localSkinFriction(sol2, p, n);
NuRe12 = nusseltNumberFn(sol2, p, n);
fprintf('Cf Re^{1/2} = %7.9f\n',CfRe12);
fprintf('Nu Re^{-1/2} = %7.9f\n',NuRe12);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the ODE function
function ff = OdeBVP (x, y)
    global p n
    denom = (p.omega * (p.delta * y(3)^2 - 1)) - n.phiMu;
    inertial = y(2)^2 - y(1)*y(3) + p.A*(y(2) + (x/2)*y(3));
    electromagnetic = n.phiSigma * p.M * y(2);
    rhsMomentum = -n.phiRho*( inertial-(n.phiBeta * p.lambda * y(4) * cos(p.alpha)));

    thermalNumer = y(2)*y(4) - y(1)*y(5) + p.A*(2*y(4) + (x/2)*y(5));
    thermalDenom = n.phiK + (4/3)*p.Rd;

    ff = [
        y(2)
        y(3)
        (rhsMomentum - electromagnetic) / denom
        y(5)
        p.Pr * n.phiRhoCp * (thermalNumer / thermalDenom)
        ];

    %Define the boundary condition
function res = OdeBC (ya, yb)
    global p n
    % res = [
    %     ya(1)-p.S
    %     ya(2)-1-(p.Sl*ya(3))
    %     n.phiK-ya(5)+(p.Bi*(1-ya(4)))
    %     yb(2)
    %     yb(4)
    %     ];


    res = [ ya(1)-p.S
        ya(2)- 1-p.Sl*ya(3)
        %ya(4)-1
        ya(5)+p.Bi*(1-ya(4))
        yb(2)
        yb(4)];

    %Setting the initial guess for first solution
function guess = OdeInit1 (x, p)
    global p
    guess =[0
        0
        0
        0
        0];

    %Setting the initial guess for second solution
function guess = OdeInit2 (x, p)

    global p

    % Softer branch-bias guess to avoid mesh blowup
    % (keeps theta' negative at eta=0 while reducing stiffness)
    guess =[p.S + 0.4*(1-exp(-x))                   % f
        1 - 1.8*exp(-x)                             % f'
        1.8*exp(-x)                                 % f''
        1 + 0.2*exp(-x)                             % theta
        -0.2*exp(-x)];                              % theta'



function CfRe12 = localSkinFriction(sol, p, n)
    y0  = deval(sol, 0);
    fpp0  = y0(3);
    CfRe12 = (n.phiMu + p.omega)*fpp0 - (1/3)*p.delta*p.omega*(fpp0^3);

function NuRe12 = nusseltNumberFn(sol, p, n)
    y0  = deval(sol, 0);
    tp0 = y0(5);
    NuRe12 = -n.phiK*(1+ (4/3)*p.Rd)*tp0;



