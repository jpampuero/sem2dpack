clearvars;
close all;

addpath D:\GitHub\sem2dpack\POST\matlab  % check this path in your SEM2DPACK installation directory

%% Plot results for a 2.5D plastic simulation
% by Ekaterina Bolotskaya
% 01/2025

%% Results directory
dir = 'psi_45_S_0.56_CF_0.63_W_10';

W       = 10e3;            % full stripe width
psi     = 45;              % stress angles
S       = 0.56;            % prestress factor
CF      = 0.63;            % closeness to failure

%% Common parameters
cp      = 5770;            % p-wave velocity
cs      = 3330;            % s-wave velocity
rho     = 2705;            % density
Dc      = 2;  
mus     = 0.6;
mud     = 0.1;
szz     = 50e6;            % normal stress
t_start = 23;              % start time for stable rupture (empirical)
t_fin   = 35;              % finish time for stable rupture (empirical)

% Compute Lame parameters
mu      = rho*cs*cs;
lambda  = rho*cp*cp-2*mu;
nu      = lambda/2/(lambda+mu);

cR      = cs*(0.862+1.14*nu)/(1+nu);      % Rayleigh wave speed

% Plotting parameters
fs      = 14;               % font size
lw      = 1.2;              % line width
isnap   = 3;                % snapshot # to plot

%% Read grid
g            = sem2d_read_specgrid(char(dir));
k            = find(g.coord(:,1)==0);
z_axis       = g.coord(k,2);
[~,idx_all]  = ismember(k,g.ibool(:));
[z_a_s,sIdx] = sort(z_axis);

%% Plastic strain
fname        = 'pla';
pla          = sem2d_snapshot_read(fname, isnap, char(dir));

% syms e11 e22 e33 e12 e13 e23
% tr  = (e11+e22+e33)/3;
% Y2  = simplify(0.5*((e11-tr)^2+(e22-tr)^2+(e33-tr)^2+2*e12^2+2*e23^2+2*e13^2));
% epl = simplify((4/3*Y2)^0.5);
    
e11   = pla.ep11;
e22   = pla.ep22;
e12   = pla.ep12;
e33   = 0;
e23   = 0;
e13   = 0;
PEEQ  = (2*(e11.^2 - e11.*e22 - e11.*e33 + 3*e12.^2 + 3*e13.^2 + e22.^2 - e22.*e33 + 3*e23.^2 + e33.^2).^(1/2))/3;
PEEQ1 = (sqrt(2)*((e11-e22).^2 + (e22-e33).^2 - (e33-e11).^2 + 6*e12.^2 + 6*e13.^2 + 6*e23.^2).^(1/2))/3;

% Plot snapshot pla (PEEQ)
figure()
sem2d_snapshot_plot(PEEQ/(szz*mus/2/mu), g);
box on;
set(gca, 'fontsize', fs-2);
xlabel('x,\ m', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('z,\ m', 'Interpreter', 'latex', 'FontSize', fs);
title('$\epsilon^{pl}_{eq}/(\tau_{s}/2\mu)$', 'Interpreter', 'latex', 'FontSize', fs+2);
caxis([0 1]);
colorbar;
    
%% Fault
data = sem2d_read_fault(strcat(char(dir),'\Flt05'));
t_or = data.dt*linspace(0,data.nt-1,data.nt);

% Plot fault
figure()
imagesc(t_or, data.x, data.st0(1)+data.st);
set(gca, 'ydir', 'normal', 'fontsize', fs-2);
xlabel('Time,\ s', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('x,\ m', 'Interpreter', 'latex', 'FontSize', fs);
title('Shear\ stress,\ Pa', 'Interpreter', 'latex', 'FontSize', fs+2);
caxis([0, szz*mus*1.1]);
colorbar;

figure()
imagesc(t_or, data.x, data.v);
hold on;
t_w = linspace(0,(data.nt-2)*data.dt, 10);
plot(t_w, t_w*cs, 'r');
plot(t_w, t_w*cR, 'r');
set(gca, 'ydir', 'normal', 'fontsize', fs-2);
xlabel('Time,\ s', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('x,\ m', 'Interpreter', 'latex', 'FontSize', fs);
title('Slip\ rate,\ m/s', 'Interpreter', 'latex', 'FontSize', fs+2);
colorbar;

figure()
imagesc(t_or, data.x, data.d);
set(gca, 'ydir', 'normal', 'fontsize', fs-2);
xlabel('Time,\ s', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('x,\ m', 'Interpreter', 'latex', 'FontSize', fs);
title('Slip,\ m', 'Interpreter', 'latex', 'FontSize', fs+2);
colorbar;

%% Rupture fronts
front_mu      = -(data.st0(1)+data.st)./(data.sn0(1)+data.sn);
front_mustick = -(data.st0(1)+data.sts)./(data.sn0(1)+data.sn);
[Trup,Tpz]    = plot_fronts_Tstick(front_mu, front_mustick, mus ,data.d, Dc, data.x, data.dt);

% Plot rupture fronts
figure()    
plot(Tpz, data.x, Trup, data.x);
legend('Tail', 'Front', 'Interpreter', 'latex', 'FontSize', fs-2, 'location', 'northwest');
title('Rupture\ fronts', 'Interpreter', 'latex', 'FontSize', fs+2);
ylabel('x,\ m', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('Time,\ s', 'Interpreter', 'latex', 'FontSize', fs);
grid on;

%% Rupture velocity
% Central difference
x_in          = linspace(min(data.x), max(data.x), fix(length(data.x)/1));
t_in          = interp1(data.x, Trup, x_in);
vrup          = zeros(size(x_in));

vrup(1)       = (x_in(2)-x_in(1))./(t_in(2)-t_in(1));
vrup(2:end-1) = (x_in(3:end)-x_in(1:end-2))./(t_in(3:end)-t_in(1:end-2));
vrup(end)     = (x_in(end)-x_in(end-1))./(t_in(end)-t_in(end-1));

vrup(vrup<0)  = 0;
vrup(vrup>cp) = cp;

frv           = mean(vrup((t_in>t_start) & (t_in<t_fin)));                 % final rupt vel
irs           = find(vrup >= frv, 1, 'first');                             % risetime

xrup          = interp1(Trup, data.x, linspace(t_start, t_fin, 100));
xpz           = interp1(Tpz, data.x, linspace(t_start, t_fin, 100));
fw            = mean(xrup-xpz);                                            % front width

fprintf('v_r = %f m/s, Risetime = %f s, Front width = %f m\n', frv, t_in(irs), fw);

% Plot rupture velocity
figure()
yline(cR);
hold on;
plot(t_in,vrup);
plot(t_in,smoothdata(vrup,'gaussian',20), 'LineWidth', lw);
yline(frv, 'g');
legend('Rayleigh speed','$v_{r}$','Smooth $v_{r}$','Final $v_{r}$', 'Interpreter', 'latex', 'FontSize', fs-2, 'location', 'southeast');
title('Rupture\ speed', 'Interpreter', 'latex', 'FontSize', fs+2);
ylabel('$v_{r}$,\ m/s', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('Time,\ s', 'Interpreter', 'latex', 'FontSize', fs);
xlim([0, max(t_w)]);
grid on;
box on;

%% Energies
Gc          = abs(szz)*(mus-mud)*Dc/2;
G0         = (abs(data.st0(1))-abs(szz*mud))^2/(mu*(pi*(1-nu)/W));

% Gc from simulation: Integrate (tau-tau_d)*slip (every node)
Gc_sim_m   = cumsum(data.v.*data.dt.*(data.st0(1)+data.st-mud*szz),2);
Gc_sim     = sum(Gc_sim_m.*data.wei(:,1),1);                             % sum over the fault length

% Heat: Integrate tau_d*slip (every node)
Heat_sim_m = cumsum(data.v.*data.dt.*(mud*szz),2);
Heat_sim   = sum(Heat_sim_m.*data.wei(:,1),1);                           % sum over the fault length

energies   = sem2d_read_and_compute_energies(dir, data.x, Trup, Gc_sim(1:end-1)', Heat_sim(1:end-1)');

%% Plot energies vs slip
% Mean slip
sl_me             = data.d;
sl_me(sl_me == 0) = NaN;
d_me              = mean(sl_me, 1, 'omitnan');
d_me_en           = interp1(t_or, d_me, energies.t_en);

figure()
yline(G0-Gc, 'k', 'LineWidth', lw*1.5);
yline(G0, 'g', 'LineWidth', lw*1.5);
yline(Gc, 'r', 'LineWidth', lw*1.5);
hold on;
plot(d_me_en, smoothdata(energies.dEp,'gaussian',20));
plot(d_me_en, smoothdata(energies.dEe,'gaussian',20));
plot(d_me_en, smoothdata(energies.dEk,'gaussian',20));
plot(d_me_en, smoothdata(energies.dHeat_sim,'gaussian',20));
plot(d_me_en, smoothdata(energies.dGc_sim,'gaussian',20));
plot(d_me_en, smoothdata(energies.dGc_sim+energies.dHeat_sim+energies.dEk+energies.dEp,'gaussian',20));
text(0.1,0.8, strcat('$\psi$ = ', num2str(psi, '%.0f'), ', S = ', num2str(S, '%.2f'), ', CF = ', num2str(CF, '%.2f'), ', W = ', num2str(W/1e3, '%.0f'), ' km'),'Units','normalized', 'Interpreter', 'latex', 'FontSize', fs);
legend('$G_0-G_c$','$G_0$','$G_c$','$dE_{pl}$','$dE_{e}$','$dE_{k}$','Heat','$G_c$ simulation','Sum = Gc+Heat+Ek+Ep',  'Interpreter', 'latex', 'FontSize', fs-2, 'location', 'southeast');
ylabel('$dE$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('Slip,\ m', 'Interpreter', 'latex', 'FontSize', fs);
axis([0 5.5 0 0.7e8]);
grid on;
box on;

%% Plot energies vs. time
figure()
yline(G0-Gc, 'k', 'LineWidth', lw*1.5);
yline(G0, 'g', 'LineWidth', lw*1.5);
yline(Gc, 'r', 'LineWidth', lw*1.5);
hold on;
plot(energies.t_en, smoothdata(energies.dEp,'gaussian',20));
plot(energies.t_en, smoothdata(energies.dEe,'gaussian',20));
plot(energies.t_en, smoothdata(energies.dEk,'gaussian',20));
plot(energies.t_en, smoothdata(energies.dHeat_sim,'gaussian',20));
plot(energies.t_en, smoothdata(energies.dGc_sim,'gaussian',20));
plot(energies.t_en, smoothdata(energies.dGc_sim+energies.dHeat_sim+energies.dEk+energies.dEp,'gaussian',20));
text(0.1,0.8, strcat('$\psi$ = ', num2str(psi, '%.0f'), ', S = ', num2str(S, '%.2f'), ', CF = ', num2str(CF, '%.2f'), ', W = ', num2str(W/1e3, '%.0f'), ' km'),'Units','normalized', 'Interpreter', 'latex', 'FontSize', fs);
legend('$G_0-G_c$','$G_0$','$G_c$','$dE_{pl}$','$dE_{e}$','$dE_{k}$','Heat','$G_c$ simulation','Sum = Gc+Heat+Ek+Ep',  'Interpreter', 'latex', 'FontSize', fs-2, 'location', 'southeast');
ylabel('$dE$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('Time,\ s', 'Interpreter', 'latex', 'FontSize', fs);
axis([0 60 0 0.7e8]);
grid on;
box on;
