% assuming TCs are on material boundaries

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
TCpos = [0, .04, .048, .056, .064, .077];

% thermal conductivities, W m-1 K-1
k_Cu = 398;
func_k_Inco = @(T) 11.45 + 1.156e-2*T + 7.72e-6*T.^2;
k_IncoUp = func_k_Inco((dat.T_Cu2 + dat.T_Inco1)/2);
k_IncoDown = func_k_Inco((dat.T_Inco2 + dat.T_Cu3)/2);
func_k_H25 = @(T) 9.9905 + .0205*T -3e-6*T.^2;
k_H25 = func_k_H25((dat.T_Inco1 + dat.T_Inco2));

% delta T
dT_IncoUp = dat.T_Cu2 - dat.T_Inco1;
dT_IncoDown = dat.T_Inco2 - dat.T_Cu3;
dT_H25 = dat.T_Inco1 - dat.T_Inco2;

% heat flux
Q_IncoUp = k_IncoUp.*dT_IncoUp./(TCpos(3)-TCpos(2));
Q_IncoDown = k_IncoDown.*dT_IncoDown./(TCpos(5)-TCpos(4));
Q_H25 = (Q_IncoUp + Q_IncoDown)/2;

% plots
figure()
hold on
grid minor
plot(dat.time, dat.T_oven)
plot(dat.time, dat.T_room)
plot(dat.time, dat.T_Cu1)
plot(dat.time, dat.T_Cu2)
plot(dat.time, dat.T_Inco1)
plot(dat.time, dat.T_Inco2)
plot(dat.time, dat.T_Cu3)
plot(dat.time, dat.T_Cu4)
xlabel('Time, s')
ylabel('Temperature, {^\circ}C')
ylim([0 120])

figure()
hold on
grid minor
plot(dat.time, Q_IncoDown)
plot(dat.time, Q_IncoUp)
xlabel('Time, s')
ylabel('Q, W m^{-2}')