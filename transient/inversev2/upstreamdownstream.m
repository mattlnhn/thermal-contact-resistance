%clear; close all; clc;

filename = 'data-16-08-2023-12-00.csv';

% geometry/material properties can be copy/pasted from main.m.

%% section geometry
% upstream copper/inco
gu.dx = .25e-3;
gu.L1 = 1.5e-3;
gu.L2 = 6.5e-3;
gu.N = (gu.L1+gu.L2)/gu.dx;
fu = ceil(gu.N*gu.L1/(gu.L1+gu.L2)):gu.N; % inco fraction
% central inco/h25/inco
gc.dx = 1.0e-4;
gc.L1 = 1.5e-3;
gc.L2 = 5.0e-3;
gc.L3 = 1.5e-3;
gc.N = (gc.L1+gc.L2+gc.L3)/gc.dx;
% downstream inco/copper
gd.dx = .25e-3;
gd.L1 = 6.5e-3;
gd.L2 = 1.5e-3;
gd.N = (gd.L1+gd.L2)/gd.dx;
fd = 1:floor(gd.N*gd.L1/(gd.L1+gd.L2)); % inco fraction

%% section material properties
% universal
Fk_Cu = @(T) 385; % [W m-1 K-1]
Fc_p_Cu = @(T) 399.5814286 - 0.0585714*T; % [J kg-1 K-1]
Frho_Cu = @(T) 8940; % [kg m-3]
Fk_Inco = @(T) 9.5164989 + 0.0216787*T + -0.0000039*T^2; % [W m-1 K-1]
Fc_p_Inco = @(T) 363.8195515 + 0.1233661*T + 0.0000527*T^2; % [J kg-1 K-1]
Frho_Inco = @(T) 8190; % [kg m-3]
Fk_H25 = @(T) 9.9905357 + 0.0205437*T - 0.000003*T^2; % [W m-1 K-1]
Fc_p_H25 = @(T) 396.5228931 + 0.2075422*T + 0.0000134*T^2; % [J kg-1 K-1]
Frho_H25 = @(T) 9070; % [kg m-3]

% upstream
mu.Fk1 = Fk_Cu;
mu.Fc_p1 = Fc_p_Cu;
mu.Frho1 = Frho_Cu;
mu.Fk2 = Fk_Inco;
mu.Fc_p2 = Fc_p_Inco;
mu.Frho2 = Frho_Inco;

% central
mc.Fk1 = Fk_Inco;
mc.Fc_p1 = Fc_p_Inco;
mc.Frho1 = Frho_Inco;
mc.Fk2 = Fk_H25;
mc.Fc_p2 = Fc_p_H25;
mc.Frho2 = Frho_H25;
mc.Fk3 = Fk_Inco;
mc.Fc_p3 = Fc_p_Inco;
mc.Frho3 = Frho_Inco;

% downstream
md.Fk1 = Fk_Inco;
md.Fc_p1 = Fc_p_Inco;
md.Frho1 = Frho_Inco;
md.Fk2 = Fk_Cu;
md.Fc_p2 = Fc_p_Cu;
md.Frho2 = Frho_Cu;

%% upstream/downstream heat flux
dt_ud = .25e-3;
dat = importInterpData(filename, dt_ud, 'makima');
timeSteps = length(dat.time);
qu = zeros(timeSteps, 1);
qd = zeros(timeSteps, 1);

% initial temp distribution
% linear interpolation between initial end temps
Tu = interp1([1; gu.N], [dat.T_Cu2(1); dat.T_Inco1(1)], 1:gu.N)';
Td = interp1([1; gd.N], [dat.T_Inco2(1); dat.T_Cu3(1)], 1:gd.N)';

for m = 1:timeSteps
    % upstream copper/inco
    Tu = temp1I_PC_TBC(dat.T_Cu2(m), dat.T_Inco1(m), Tu, dt_ud, gu, mu);
    coefu = [ones(length(fu), 1) gu.dx*fu']\Tu(fu);
    qu(m) = -mu.Fk2(.5*(dat.T_Cu2(m)+dat.T_Inco1(m)))*coefu(2);

    % downstream inco/copper
    Td = temp1I_PC_TBC(dat.T_Inco2(m), dat.T_Cu3(m), Td, dt_ud, gd, md);
    coefd = [ones(length(fd), 1) gd.dx*fd']\Td(fd);
    qd(m) = -md.Fk1(.5*(dat.T_Inco2(m)+dat.T_Cu3(m)))*coefd(2);

    fprintf('%.2f%% complete.\n', 100*m/timeSteps)
end

save upstreamdownstream.mat qu qd dt_ud