clear; close all; clc;

filename = "5N.dat";

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

%% import
dat = readtable(filename);
totalSteps = length(dat.time);
qu = zeros(totalSteps+1, 2);
qd = qu;

%% initial temp distribution
% linear interpolation between initial end temps
initialT = interp1([1; gu.N], [dat.T_Cu2(1); dat.T_Inco1(1)], 1:gu.N)';
%Td = interp1([1; gd.N], [dat.T_Inco2(1); dat.T_Cu3(1)], 1:gd.N)';

%% parameters
dt = .25e-3;    % finite volume time step
epsilon = 1e-2; % small fraction for newton method
TOLdq = 1e-3;   % tolerance for relative change in dq
TOLerror = 1e-9;    % tolerance for relative change in error
maxIter = 25;   % max iterations
qInitial = 1000;    % first guess at q
qu(1, :) = qInitial;
qd(1, :) = qInitial;

%% time iteration upstream
for m = 1:totalSteps-1
    % reset counters
    converged = 0;
    iter = 1;
    errorL = 1e20;
    errorR = errorL;
    dqL = errorL;
    dqR = errorL;

    currentTime = dat.time(m);
    nextTime = dat.time(m+1);
    currentSteps = round((nextTime-currentTime)/dt);

    T = zeros(gu.N, currentSteps+1);    % temp matrices incl. initial time
    T(:, 1) = initialT;
    TdqL = T;
    TdqR = T;

    qL = qu(m, 1);
    qR = qu(m, 2);

    Y = [dat.T_Cu2(m);
        dat.T_Inco1(m)];
    Tavg = .5*(Y(1) + Y(2));

    while converged == 0
        prevdqL = dqL;
        prevdqR = dqR;
        preverrorL = errorL;
        preverrorR = errorR;

        for n = 1:currentSteps
            T(:, n+1) = temp1I_PC_QBC(qL, qR, Tavg, T(:, n), dt, gu, mu);
            TdqL(:, n+1) = temp1I_PC_QBC((1+epsilon)*qL, qR, Tavg, ...
                TdqL(:, n), dt, gu, mu);
            TdqR(:, n+1) = temp1I_PC_QBC(qL, (1+epsilon)*qR, Tavg, ...
                TdqR(:, n), dt, gu, mu);
        end

        % T at sensor location, excl. initial
        Ts = [T(1, 2:end); T(end, 2:end)];
        TdqLs = [TdqL(1, 2:end); TdqL(end, 2:end)];
        TdqRs = [TdqR(1, 2:end); TdqR(end, 2:end)];

        % sensitivity coefficients
        XL = (TdqLs - Ts)./(epsilon*qL);
        XR = (TdqRs - Ts)./(epsilon*qR);

        % newton method increment
        dqL = sum((Y-Ts).*XL, 'all')/sum(XL.^2, 'all');
        dqR = sum((Y-Ts).*XR, 'all')/sum(XR.^2, 'all');

        qL = qL + dqL;
        qR = qR + dqR;

        errorL = sum(.5*(Y-TdqLs).^2, 'all');
        errorR = sum(.5*(Y-TdqRs).^2, 'all');

        if abs((dqL-prevdqL)/prevdqL) < TOLdq && ...
                abs((dqR-prevdqR)/prevdqR) < TOLdq
            converged = 1;
            fprintf('TOLdq. ')
        end
        if abs((errorL-preverrorL)/preverrorL) < TOLerror && ...
                abs((errorR-preverrorR)/preverrorR) < TOLerror
            converged = 1;
            fprintf('TOLerror. ')
        end
        if iter > maxIter
            converged = 1;
            fprintf('Max iter. ')
        end

        iter = iter + 1;

    end

    % save q
    qu(m+1, :) = [qL qR];
    
    % next step initial temp
    for n = 1:currentSteps
        initialT = temp1I_PC_QBC(qL, qR, Tavg, initialT, dt, gu, mu);
    end

    fprintf('%.3f%% complete.\n', 100*m/totalSteps)
end

qu = qu(:, 2);
qd = qd(:, 1);
newfilename = filename+".mat";
save(newfilename, "qu", "qd", "dt_ud")

