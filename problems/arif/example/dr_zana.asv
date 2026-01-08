function CodeA

format long g 

%Define all parameters 
global phi1 phi2 phi Om Del Pr A L alpha Rd Bi Ec lambda A1 A2 A3 A4 A5 A6 M S

%Boundary layer thickness & stepsize
etaMin = 0;
etaMax1 = 20; %blt 1st solution 20
etaMax2 = 20; %blt 2nd solution 20
stepsize1 = etaMax1;
stepsize2= etaMax2;

% Input for the parameters
phi1 = 0.01; phi2 = 0.01;  % Nanoparticle volume fraction
Pr = 6.2;  % Prandtl number Engine oil
lambda = -1; % mixed parameter
Om = 0.1; Del = 0.2; % Eyring powell fluid parameter
L = 0.2;  % Slip parameter
M = 0.2; % Magnetic parameter
S = 2; % Suction parameter
A = -1.2; %unsteady parameter
alpha = pi/4; %angle
Ec = 0.2; % eckert
Rd = 0.2; % radiation
Bi = 0.2; %Biot number

phi = phi1+phi2;

rho_1 = 3970; C_1 = 765; k_1 = 40; s_1 = 3.5*10^7; b_1 = 0.85*(10^-5); %Alumina
rho_2 = 8933; C_2 = 385; k_2 = 401; s_2 = 5.96*10^7; b_2 = 1.67*(10^-5); %Copper
rho_f = 997.1; C_f = 4179; k_f = 0.613; s_f = 0.05; b_f = 21; %water

A1 = 1/((1-phi)^2.5);
A2 = (1-phi)+phi1*(rho_1/rho_f)+phi2*(rho_2/rho_f);
A3 = (1-phi)+phi1*(b_1/b_f)+phi2*(b_2/b_f);
A4 = (((phi1*s_1+phi2*s_2)/phi)+2*s_f+2*(phi1*s_1+phi2*s_2)-2*phi*s_f)/(((phi1*s_1+phi2*s_2)/phi)+2*s_f-(phi1*s_1+phi2*s_2)+phi*s_f);
A5 = (1-phi)+phi1*(rho_1*C_1)/(rho_f*C_f)+phi2*(rho_2*C_2)/(rho_f*C_f);
A6 = (((phi1*k_1+phi2*k_2)/phi)+2*k_f+2*(phi1*k_1+phi2*k_2)-2*phi*k_f)/(((phi1*k_1+phi2*k_2)/phi)+2*k_f-(phi1*k_1+phi2*k_2)+phi*k_f);

%%%%%%%%%%%%%%%%%%%%%%   first solution   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    options = bvpset('stats','off','RelTol',1e-10); 
    solinit = bvpinit (linspace (etaMin, etaMax1, stepsize1), @OdeInit1);
    sol = bvp4c (@OdeBVP, @OdeBC, solinit, options);
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
    fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%   second solution   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    options = bvpset('stats','off','RelTol',1e-10); 
    solinit = bvpinit (linspace (etaMin, etaMax2, stepsize2), @OdeInit2);
    sol = bvp4c (@OdeBVP, @OdeBC, solinit, options);
    eta = linspace (etaMin, etaMax2, stepsize2);
    y = deval (sol, eta);
     
    figure(1)
    plot(sol.x,sol.y(2,:),'--r')
    xlabel('\eta')
    ylabel('f`(\eta)')
    hold on
    
    figure(2)
    plot(sol.x,sol.y(4,:),'--r')
    xlabel('\eta')
    ylabel('t(\eta)')
    hold on
    
    %Saving the output in txt file for second solution
    descris = [sol.x; sol.y];
    save 'lower.txt' descris -ascii
     
    %Displaying the output for second solution 
    fprintf('\nSecond solution:\n');
    fprintf('f"(0) = %7.9f\n',y(3));
    fprintf('-t`(0) = %7.9f\n',-y(5));
    fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%Define the ODE function
function ff = OdeBVP (x, y, Pr, A1, A2, A3, A4, A5, A6, Om, Del, M, lambda, A, alpha, Ec, Rd)

global Pr A1 A2 A3 A4 A5 A6 Om Del M lambda A alpha Ec Rd

 ff = [y(2)   
       y(3)   
       (1/(Om*(Del*y(3)*y(3)-1)-A1))*(-A2*(y(2)*y(2)-y(1)*y(3)+A*(y(2)+(x/2)*y(3))-A3*lambda*y(4)*cos(alpha))-A4*M*y(2))
       y(5)   
       Pr*A5*(1/(A6+(4/3)*Rd))*(-(A1/A5)*Ec*y(3)*y(3)+y(2)*y(4)-y(1)*y(5)+A*(2*y(4)+(x/2)*y(5)))]; 

%Define the boundary condition
function res = OdeBC (ya, yb, L, S, Bi)

global L S Bi
res = [ ya(1)-S
        ya(2)-1-L*ya(3)
        %ya(4)-1
        ya(5)+Bi*(1-ya(4))
        yb(2)
        yb(4)];           
 
%Setting the initial guess for first solution  
function v = OdeInit1 (x, L)

global L
    v =[0               
     0                     
     0                         
     0               
     0];
 
%Setting the initial guess for second solution 
function v1 = OdeInit2 (x, L)

global L

    v1 =[0.2+exp(-x)                  % f
     exp(-x)                        % f'
     -exp(-x)                         % f''
     exp(-x)               % θ  ≈1 at η=0
     exp(-x)];
