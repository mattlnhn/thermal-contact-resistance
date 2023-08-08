% material boundaries in correct locations

clear; close all;

% temp data
dat = readtable('tcc_data.csv');
% dat.time      = time
% dat.T_oven    = oven temperature
% dat.T_Cu1     = TC1, x = 0 mm
% dat.T_Cu2     = TC2, x = 40 mm
% dat.T_Inco1   = TC3, x = 48 mm
% dat.T_Inco2   = TC4, x = 56 mm
% dat.T_Cu3     = TC5, x = 64 mm
% dat.T_Cu4     = TC6, x = 77 mm
% dat.T_room    = room temperature
TCpos = [0, .04, .048, .056, .064, .077]; % m

% thermal conductivities, W m-1 K-1
k_Cu = 398; % W m-2 K-1
func_k_Inco = @(T) 11.45 + 1.156e-2*T + 7.72e-6*T.^2;
k_IncoUp = func_k_Inco((dat.T_Cu2 + dat.T_Inco1)/2); % W m-2 K-1
k_IncoDown = func_k_Inco((dat.T_Inco2 + dat.T_Cu3)/2); % W m-2 K-1
func_k_H25 = @(T) 9.9905 + .0205*T -3e-6*T.^2;
k_H25 = func_k_H25((dat.T_Inco1 + dat.T_Inco2)/2); % W m-2 K-1

% delta T
dT_IncoUp = dat.T_Cu2 - dat.T_Inco1; % K
dT_IncoDown = dat.T_Inco2 - dat.T_Cu3; % K
dT_H25 = dat.T_Inco1 - dat.T_Inco2; % K

% heat fluxes
t_Cu = 1.5e-3; % thickness, m
t_Inco = 6.5e-3; % m
h_IncoUp = ((t_Cu/k_Cu) + (t_Inco./k_IncoUp)).^-1; % htc, W m-2 K-1
h_IncoDown = ((t_Cu/k_Cu) + (t_Inco./k_IncoDown)).^-1; % W m-2 K-1
Q_IncoUp = h_IncoUp.*dT_IncoUp; % heat flux, W m-2
Q_IncoDown = h_IncoDown.*dT_IncoDown; % W m-2
Q_H25 = (Q_IncoUp + Q_IncoDown)/2; % W m-2

% resistance R = 1/h
t_Inco = 1.5e-3; % m
t_H25 = 5.0e-3; % m
h_cond = (t_Inco./k_IncoUp) + (t_H25./k_H25) + (t_Inco./k_IncoDown); % W m-2 K-1
R = .5*((dT_H25./Q_H25)-h_cond); % m2 K W-1
h = R.^-1; % W m-2 K-1

% plots
figure()
hold on
grid minor
plot(dat.time, dat.T_oven)
plot(dat.time, dat.T_Cu1)
plot(dat.time, dat.T_Cu2)
plot(dat.time, dat.T_Inco1)
plot(dat.time, dat.T_Inco2)
plot(dat.time, dat.T_Cu3)
plot(dat.time, dat.T_Cu4)
plot(dat.time, dat.T_room)
xlabel('Time, s')
ylabel('Temperature, {^\circ}C')
ylim([0 120])
legend({'oven', 'Cu1', 'Cu2', 'Inco1', 'Inco2', 'Cu3', 'Cu4', 'room'})
cmap = turbo(8); % colours
ax = gca;
ax.ColorOrder = cmap;

figure()
hold on
grid minor
plot(dat.time, Q_IncoDown, 'r')
plot(dat.time, Q_IncoUp, 'b')
xlabel('Time, s')
ylabel('Q, W m^{-2}')
legend({'Q_{IncoDown}', 'Q_{IncoUp}'})

figure()
hold on
grid minor
plot(dat.time, R, 'r')
xlabel('Time, s')
ylabel('Contact resistance, m^2 K W^{-1}')
ylim([0 1.5e-3])

