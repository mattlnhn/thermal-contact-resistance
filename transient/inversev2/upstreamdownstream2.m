clear; close all; clc;

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
dat = importInterpData(filename, dt_ud, 'makima', 'gaussian', 1000);
timeSteps = length(dat.time);
qu = zeros(timeSteps, 2);
qd = zeros(timeSteps, 2);

% initial temp distribution
% linear interpolation between initial end temps
Tu = interp1([1; gu.N], [dat.T_Cu2(1); dat.T_Inco1(1)], 1:gu.N)';
Td = interp1([1; gd.N], [dat.T_Inco2(1); dat.T_Cu3(1)], 1:gd.N)';

r = 5;
epsilon = 1e-6;
TOLdq = 1e-3;
TOLerror = 1e-3;
qInitial = 1000;

Tcurrent = zeros(gu.N, r);
TqLcurrent = zeros(gu.N, r);
TdqRcurent = zeros(gu.N, r);
T = zeros(2, r);
TdH = zeros(2, r);
X = zeros(2, r);

for m = 1:timeSteps
    %% upstream copper/inco
    if m ~= 1
        qL = qu(m-1, 1);
        qR = qu(m-1, 2);
    else
        qL = qInitial;
        qR = qInitial;
    end

    converged = 0;
    iter = 1;
    errorL = 1e6;
    errorR = 1e6;
    dqL = 1e6;
    dqR = 1e6;

    Yu = [dat.T_Cu2(m:m+r-1); dat.T_Inco1(m:m+r-1)];
    Tavg = .5*(Yu(1, :) + Yu(2, :));

    while converged == 0
        prevdqL = dqL;
        prevdqR = dqR;
        preverrorL = errorL;
        preverrorR = errorR;

        eqL = epsilon*qL;
        eqR = epsilon*qR;

        Tcurrent(:, 1) = temp1I_PC_QBC(qL, qR, Tavg(1), Tu, dt_ud, gu, mu);
        TdqLcurrent(:, 1) = temp1I_PC_QBC(qL+eqL, qR, Tavg(1), Tu, dt_ud, gu, mu);
        TdqRcurrent(:, 1) = temp1I_PC_QBC(qL, qR+eqR, Tavg(1), Tu, dt_ud, gu, mu);

        for j = 2:r
            Tcurrent(:, j) = temp1I_PC_QBC(qL, qR, Tavg(j), Tcurrent(:, j-1), dt_ud, gu, mu);
            TdqLcurrent(:, j) = temp1I_PC_QBC(qL+eqL, qR, Tavg(j), TdqLcurrent(:, j-1), dt_ud, gu, mu);
            TdqRcurrent(:, j) = temp1I_PC_QBC(qL, qR+eqR, Tavg(j), TdqRcurrent(:, j-1), dt_ud, gu, mu);
        end

        T = [Tcurrent(1, :); Tcurrent(end, :)];
        TdqL = [TdqLcurrent(1, :); TdqLcurrent(end, :)];
        TdqR = [TdqRcurrent(1, :); TdqRcurrent(end, :)];
        XL = (TdqL - T)./eqL;
        XR = (TdqR - T)./eqR;

        dqL = sum((Yu-T).*XL, 'all')/sum(XL.^2, 'all');
        qL = qL + dqL;
        dqR = sum((Yu-T).*XR, 'all')/sum(XR.^2, 'all');
        qR = qR + dqR;

        errorL = sum(.5*(Yu-T-dqL.*XL).^2, 'all');
        errorR = sum(.5*(Yu-T-dqR.*XR).^2, 'all');

        if abs((dqL-prevdqL)/prevdqL) < TOLdq && abs((dqR-prevdqR)/prevdqR) < TOLdq
            converged = 1;
            fprintf('U: TOLdq. ')
        end
        if abs((errorL-preverrorL)/preverrorL) < TOLerror && abs((errorR-preverrorR)/preverrorR) < TOLerror
            converged = 1;
            fprintf('U: TOLerror. ')
        end
        if iter > 25
            converged = 1;
            fprintf('U: Max iter. ')
        end

        % to deal with startup error
        % in this set up we know -ve heat flux is unphysical
        if qL < 0
            qL = qInitial;
            fprintf('Caught qL < 0. ')
        end
        if qR < 0
            qR = qInitial;
            fprintf('Caught qR < 0. ')
        end

        Tu = Tcurrent(:, 1);
        qu(m, :) = [qL qR];

        iter = iter + 1;
    end

    %% downstream inco/copper
    if m ~= 1
        qL = qd(m-1, 1);
        qR = qd(m-1, 2);
    else
        qL = qInitial;
        qR = qInitial;
    end

    converged = 0;
    iter = 1;
    errorL = 1e6;
    errorR = 1e6;
    dqL = 1e6;
    dqR = 1e6;

    Yu = [dat.T_Inco2(m:m+r-1); dat.T_Cu3(m:m+r-1)];
    Tavg = .5*(Yu(1, :) + Yu(2, :));

    while converged == 0
        prevdqL = dqL;
        prevdqR = dqR;
        preverrorL = errorL;
        preverrorR = errorR;

        eqL = epsilon*qL;
        eqR = epsilon*qR;

        Tcurrent(:, 1) = temp1I_PC_QBC(qL, qR, Tavg(1), Td, dt_ud, gu, mu);
        TdqLcurrent(:, 1) = temp1I_PC_QBC(qL+eqL, qR, Tavg(1), Td, dt_ud, gu, mu);
        TdqRcurrent(:, 1) = temp1I_PC_QBC(qL, qR+eqR, Tavg(1), Td, dt_ud, gu, mu);

        for j = 2:r
            Tcurrent(:, j) = temp1I_PC_QBC(qL, qR, Tavg(j), Tcurrent(:, j-1), dt_ud, gu, mu);
            TdqLcurrent(:, j) = temp1I_PC_QBC(qL+eqL, qR, Tavg(j), TdqLcurrent(:, j-1), dt_ud, gu, mu);
            TdqRcurrent(:, j) = temp1I_PC_QBC(qL, qR+eqR, Tavg(j), TdqRcurrent(:, j-1), dt_ud, gu, mu);
        end

        T = [Tcurrent(1, :); Tcurrent(end, :)];
        TdqL = [TdqLcurrent(1, :); TdqLcurrent(end, :)];
        TdqR = [TdqRcurrent(1, :); TdqRcurrent(end, :)];
        XL = (TdqL - T)./eqL;
        XR = (TdqR - T)./eqR;

        dqL = sum((Yu-T).*XL, 'all')/sum(XL.^2, 'all');
        qL = qL + dqL;
        dqR = sum((Yu-T).*XR, 'all')/sum(XR.^2, 'all');
        qR = qR + dqR;

        errorL = sum(.5*(Yu-T-dqL.*XL).^2, 'all');
        errorR = sum(.5*(Yu-T-dqR.*XR).^2, 'all');

        if abs((dqL-prevdqL)/prevdqL) < TOLdq && abs((dqR-prevdqR)/prevdqR) < TOLdq
            converged = 1;
            fprintf('D: TOLdq. ')
        end
        if abs((errorL-preverrorL)/preverrorL) < TOLerror && abs((errorR-preverrorR)/preverrorR) < TOLerror
            converged = 1;
            fprintf('D: TOLerror. ')
        end
        if iter > 25
            converged = 1;
            fprintf('D: Max iter. ')
        end

        % to deal with startup error
        % in this set up we know -ve heat flux is unphysical
        if qL < 0
            qL = qInitial;
            fprintf('Caught qL < 0. ')
        end
        if qR < 0
            qR = qInitial;
            fprintf('Caught qR < 0. ')
        end

        Td = Tcurrent(:, 1);
        qd(m, :) = [qL qR];

        iter = iter + 1;
    end

    fprintf('%.3f%% complete.\n', 100*m/timeSteps)
end

qu = qu(:, 2);
qd = qd(:, 1);
save upstreamdownstream2.mat qu qd dt_ud

