%% Script to compute initial strains and stresses for sem2dpack input file
%
% This script calculates the initial stress and strain in the bulk and on
% the fault based on a given set of parameters related to the fault geometry,
% material properties, and frictional behavior. 
% Computed values for initial strains, fault shear stress and cohesion 
% are then used as input values for a sem2dpack input file.
%
% Mohr circle representation to visualize the stress states is included
%    (initial and post-slip) and plasticity effects.
%
% Input:
% - Fault orientation (angle of principal stress - psi)
% - Material properties (e.g., density, wave velocities, friction coefficients)
% - Fault stress parameters (e.g., sigma_{zz}, S - prestress, CF - closeness to failure)
%
% Output:
% - Initial stress and strain states printed to the console.
% - Visualization of Mohr circles for the fault stress state.
% - Energies (G_c, G_0), and their ratio) and cohesion values.
%
% by Ekaterina Bolotskaya (modified from J-P Ampuero)
% 01/2025

clearvars;
close all;

%% Input values
psi    = 45;              % direction of max principal stress to fault (degrees)
S      = 0.56;            % strength excess parameter = (taus-tau0)/(tau0-taud)
szz    = -50e6;           % fault normal stress (sigma)
W      = 10e3;
CF     = 0.63;            % closeness to failure (Mohr circle radius / distance to yield surface)

%% Material properties
cp     = 5770;            % p-wave velocity
cs     = 3330;            % s-wave velocity
rho    = 2705;            % density
phi    = 30;              % internal friction angle (degrees)

%% Frictional properties (SW friction)
mus    = 0.6;            % static friction coefficient
mud    = 0.1;            % dynamic friction coefficient
Dc     = 2;              % slip-weakening distance

%% Plotting variables
lw = 1.2;
fs = 11;

% Compute Lame parameters
mu     = rho*cs*cs;
lambda = rho*cp*cp-2*mu;
nu     = lambda/2/(lambda+mu);
E      = mu * (3*lambda + 2*mu) / (lambda + mu);

%% Initial stress
sxz    = (mus+S*mud)/(1+S)*(-szz);
sxx    = szz-2*sxz/tand(2*psi);
syy    = nu*(sxx+szz);       % plane strain
% syy    = 0.5*(sxx+szz);    % alternative assumption
s0      = [sxx szz sxz];
fprintf('Initial stress = %e, %e, %e \n',s0);
fprintf('  Fault shear stress = %e \n',sxz);
fprintf('  Fault shear/normal stress = %e \n',sxz/szz);

%% Initial strain
% exz    = sxz/(2*mu);
% exx    = sxx/(2*mu) - lambda/(2*mu*(2*mu+3*lambda))*(sxx+szz+syy);
% ezz    = szz/(2*mu) - lambda/(2*mu*(2*mu+3*lambda))*(sxx+szz+syy);
% eyy    = syy/(2*mu) - lambda/(2*mu*(2*mu+3*lambda))*(sxx+szz+syy); % should be zero

exz    = (1+nu)/E*sxz;
exx    = (1+nu)/E*((1-nu)*sxx - nu*szz);
ezz    = (1+nu)/E*(-nu*sxx + (1-nu)*szz);
eyy    = syy/(2*mu) - lambda/(2*mu*(2*mu+3*lambda))*(sxx+szz+syy); % should be zero

e0      = [exx ezz exz];
fprintf('Initial strain (exx, ezz, exz) = %e, %e, %e \n',e0);

%% Check initial stresses from Fortran mat_plastic
s0f = zeros(3, 1); 

s0f(1) = (lambda + 2*mu) * e0(1) + lambda * e0(2); 
s0f(2) = lambda * e0(1) + (lambda + 2*mu) * e0(2);  
s0f(3) = 2*mu * e0(3);  

fprintf('Initial stress computed like in mat_plastic.f90 = %e, %e, %e \n',s0f);

%% Energies
Gc     = abs(szz)*(mus-mud)*Dc/2;
G0     = (abs(sxz)-abs(szz*mud))^2/(mu*(pi*(1-nu)/W));    % buried fault gamma = 1/pi
fprintf('G_c = %f, G_0 = %f, G_c/G_0 = %f \n', Gc, G0, Gc/G0);

%% "Closeness to failure" - max shear/distance to yield surface
syms coh

c_val  = eval(solve(sqrt(sxz.^2+(sxx-szz).^2/4)/(-(sxx+szz)/2*sind(phi)+coh*cosd(phi)) == CF, coh));
fprintf('CF = %f Cohesion = %e \n', CF, c_val);

[s_m, tau_m, cent_c] = mohr(sxx, szz, sxz, psi, 500);
% Update stresses
sxzu    = 2*(mus*(-szz)-sxz)+sxz;
sxxu    = szz-2/tand(2*psi)*sxzu;
[s_mu, tau_mu, cent_cu] = mohr(sxxu, szz, sxzu, psi, 500);
    
figure()
plot(s_m, tau_m, 'k', 'Linewidth', lw);
hold on;
plot(s_mu, tau_mu, 'k--', 'Linewidth', lw);
plot(linspace(0,min([s_mu s_m]), 100), -linspace(0,min([s_mu s_m]), 100)*tand(phi)+c_val, 'r', 'Linewidth', lw);
plot(szz, sxz, 'r*', 'Linewidth', lw);
axis equal tight;
ylim([0 -min([s_mu s_m])*tand(phi)+c_val]);
grid on;
xlabel('$\sigma_{n}$,\ Pa', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\tau$,\ Pa', 'Interpreter', 'latex', 'FontSize', fs);
legend('Initial','Final?','Plasticity', 'Fault', 'Interpreter', 'latex', 'FontSize', fs);


%% Characteristic lengths
% Lnuc = 1.158*mu/(1-nu)*Dc/(-szz*(mus-mud));
% Lpz = 9*pi/16 *mu/(1-nu) *0.5*Dc/(-szz*(mus-mud));
% fprintf('The mesh should resolve well the following lengths:\n');
% fprintf('  Lnuc = %g (nucleation length from Uenishi and Rice (2003))\n',Lnuc);
% fprintf('  Lpz  = %g (static process zone length)\n',Lpz);
% 
% Lc = 4/pi *mu/(1-nu) *0.5*Dc*(-szz)*(mus-mud) / (abs(sxz)+szz*mud)^2 ;
% fprintf('The nucleation region should be larger than:\n');
% fprintf('  Lc   = %g (critical crack length)\n',Lc);

%% Mohr circle function
function [sigma_mohr, tau_mohr, center_circle] = mohr(sigma_x, sigma_y, tau_xy, psi, gridsize)
phi=linspace(-psi/180*pi,pi/2-psi/180*pi,gridsize);
sigma_mohr=(sigma_x+sigma_y)/2+(sigma_x-sigma_y)/2*cos(2*phi)+...
    tau_xy*sin(2*phi);
tau_mohr=-(sigma_x-sigma_y)/2*sin(2*phi)+...
    tau_xy*cos(2*phi);
% sigma_1=(sigma_x+sigma_y)/2-sqrt(((sigma_x-sigma_y)/2)^2+tau_xy^2);
% sigma_2=(sigma_x+sigma_y)/2+sqrt(((sigma_x-sigma_y)/2)^2+tau_xy^2);
% tau_1=sqrt(((sigma_x-sigma_y)/2)^2+tau_xy^2);
% tau_2=-tau_1;
center_circle=(sigma_x+sigma_y)/2;
%phi_p=atan(2*tau_xy/(sigma_x-sigma_y))/2;
end
