%clear; close all; clc;

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

%% parameters
dt = .25e-3;    % finite volume time step

r = 25; % future points of data considered
epsilon = 1e-2; % small fraction for newton method
RTOLq = 1e-3;    % tolerance for relative change in q
RTOLdq = 1e-2;  % tolerance for relative change in dq
RTOLerror = 1e-6;   % tolerance for relative change in error
maxIter = 20;   % max iterations

qInitial = 10000;   % first guess at q
qu(1, :) = qInitial;
qd(1, :) = qInitial;

%% time iteration upstream
% initial temp distribution
initialT = interp1([1; gu.N], [dat.T_Cu2(1); dat.T_Inco1(1)], 1:gu.N)';

for m = 1:totalSteps-1
    % reset counters
    converged = 0;
    iter = 1;
    errorL = 1e20;
    errorR = errorL;
    dqL = errorL;
    dqR = errorL;

    % time data
    time = dat.time(m:m+r);
    steps = round((time - time(1))/dt);

    % temp matrices incl. initial time
    T = zeros(gu.N, steps(end)+1);
    T(:, 1) = initialT;
    TdqL = T;
    TdqR = T;

    qL = qu(m, 1);
    qR = qu(m, 2);

    Y = [dat.T_Cu2(m:m+r)';
        dat.T_Inco1(m:m+r)'];
    Tavg = mean(Y, 'all');

    while converged == 0
        prevdqL = dqL;
        prevdqR = dqR;
        preverrorL = errorL;
        preverrorR = errorR;

        for n = 1:steps(end)
            T(:, n+1) = temp1I_PC_QBC(qL, qR, Tavg, T(:, n), dt, gu, mu);
            TdqL(:, n+1) = temp1I_PC_QBC((1+epsilon)*qL, qR, Tavg, ...
                TdqL(:, n), dt, gu, mu);
            TdqR(:, n+1) = temp1I_PC_QBC(qL, (1+epsilon)*qR, Tavg, ...
                TdqR(:, n), dt, gu, mu);
        end

        % T at sensor locations & measurement times, excl. initial
        Ts = [T(1, steps+1); T(end, steps+1)];
        TdqLs = [TdqL(1, steps+1); TdqL(end, steps+1)];
        TdqRs = [TdqR(1, steps+1); TdqR(end, steps+1)];

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

        % convergence criteria
        % relative change in q
        if abs(dqL/qL) < RTOLq && abs(dqR/qR) < RTOLq
            converged = 1;
            fprintf('U: RTOLq. ')
        end
        % relative change in dq
        if abs((dqL-prevdqL)/prevdqL) < RTOLdq && ...
                abs((dqR-prevdqR)/prevdqR) < RTOLdq
            converged = 1;
            fprintf('U: RTOLdq. ')
        end
        % relative change in error
        if abs((errorL-preverrorL)/preverrorL) < RTOLerror && ...
                abs((errorR-preverrorR)/preverrorR) < RTOLerror
            converged = 1;
            fprintf('U: RTOLerror. ')
        end
        % timeout
        if iter > maxIter
            converged = 1;
            fprintf('U: Max iter. ')
        end

        iter = iter + 1;

    end

    % save q
    qu(m+1, :) = [qL qR];
    
    % next step initial temp
    for n = 1:steps(2)
        initialT = temp1I_PC_QBC(qL, qR, Tavg, initialT, dt, gu, mu);
    end

    fprintf('%.3f%% complete.\n', 100*m/totalSteps)
end

%% time iteration downstream
% initial temp distribution
initialT = interp1([1; gd.N], [dat.T_Inco2(1); dat.T_Cu3(1)], 1:gd.N)';

for m = 1:totalSteps-1
    % reset counters
    converged = 0;
    iter = 1;
    errorL = 1e20;
    errorR = errorL;
    dqL = errorL;
    dqR = errorL;

    % time data
    time = dat.time(m:m+r);
    steps = round((time - time(1))/dt);

    % temp matrices incl. initial time
    T = zeros(gd.N, steps(end)+1);
    T(:, 1) = initialT;
    TdqL = T;
    TdqR = T;

    qL = qd(m, 1);
    qR = qd(m, 2);

    Y = [dat.T_Inco2(m:m+r)';
        dat.T_Cu3(m:m+r)'];
    Tavg = mean(Y, 'all');

    while converged == 0
        prevdqL = dqL;
        prevdqR = dqR;
        preverrorL = errorL;
        preverrorR = errorR;

        for n = 1:steps(end)
            T(:, n+1) = temp1I_PC_QBC(qL, qR, Tavg, T(:, n), dt, gd, md);
            TdqL(:, n+1) = temp1I_PC_QBC((1+epsilon)*qL, qR, Tavg, ...
                TdqL(:, n), dt, gd, md);
            TdqR(:, n+1) = temp1I_PC_QBC(qL, (1+epsilon)*qR, Tavg, ...
                TdqR(:, n), dt, gd, md);
        end

        % T at sensor locations & measurement times, excl. initial
        Ts = [T(1, steps+1); T(end, steps+1)];
        TdqLs = [TdqL(1, steps+1); TdqL(end, steps+1)];
        TdqRs = [TdqR(1, steps+1); TdqR(end, steps+1)];

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

        % convergence criteria
        % relative change in q
        if abs(dqL/qL) < RTOLq && abs(dqR/qR) < RTOLq
            converged = 1;
            fprintf('D: RTOLq. ')
        end
        % relative change in dq
        if abs((dqL-prevdqL)/prevdqL) < RTOLdq && ...
                abs((dqR-prevdqR)/prevdqR) < RTOLdq
            converged = 1;
            fprintf('D: RTOLdq. ')
        end
        % relative change in error
        if abs((errorL-preverrorL)/preverrorL) < RTOLerror && ...
                abs((errorR-preverrorR)/preverrorR) < RTOLerror
            converged = 1;
            fprintf('D: RTOLerror. ')
        end
        % timeout
        if iter > maxIter
            converged = 1;
            fprintf('D: Max iter. ')
        end

        iter = iter + 1;

    end

    % save q
    qd(m+1, :) = [qL qR];
    
    % next step initial temp
    for n = 1:steps(2)
        initialT = temp1I_PC_QBC(qL, qR, Tavg, initialT, dt, gd, md);
    end

    fprintf('%.3f%% complete.\n', 100*m/totalSteps)
end

qu = qu(:, 2);
qd = qd(:, 1);
newfilename = filename+".mat";
save(newfilename, "qu", "qd", "dt_ud")

