%clear; close all; clc;
clear;
filename = "20230830_1_5Ncycle8.dat";


%% section geometry
% geometry/material properties can be copy/pasted from main.m.
% upstream copper/inco
% rows: dx, [L1, L2, ...], [N1, N2, ...]
gu = cell(3, 1);
gu{1} = .25e-3;
gu{2} = [1.5e-3, 6.5e-3];
gu{3} = [gu{2}(1)/gu{1}, gu{2}(2)/gu{1}];
Nu = round(sum(gu{3}, "all"));
% central inco/h25/inco
% rows: dx, [L1, L2, ...], [N1, N2, ...]
gc = cell(4, 1);
gc{1} = 1.0e-4;
gc{2} = [1.5e-3, 5.0e-3, 1.5e-3];
gc{3} = [gc{2}(1)/gc{1}, gc{2}(2)/gc{1}, gc{2}(3)/gc{1}];
Nc = round(sum(gc{3}, "all"));
TCH25 = floor(Nc*3e-3/8e-3); % node # for H25 thermocouple
% downstream inco/copper
% rows: dx, [L1, L2, ...], [N1, N2, ...]
gd = cell(3, 1);
gd{1} = .25e-3;
gd{2} = [6.5e-3, 1.5e-3];
gd{3} = [gd{2}(1)/gd{1}, gd{2}(2)/gd{1}];
Nd = round(sum(gd{3}, "all"));


%% section A matrices
Au = diag(ones(Nu-1, 1), -1) + diag(-2*ones(Nu, 1), 0) + ...
    diag(ones(Nu-1, 1), 1);
Au(1, 1) = -1; Au(end, end) = -1;
Ac = diag(ones(Nc-1, 1), -1) + diag(-2*ones(Nc, 1), 0) + ...
    diag(ones(Nc-1, 1), 1);
Ac(1, 1) = -1; Ac(end, end) = -1;
Ad = diag(ones(Nd-1, 1), -1) + diag(-2*ones(Nd, 1), 0) + ...
    diag(ones(Nd-1, 1), 1);
Ad(1, 1) = -1; Ad(end, end) = -1;


%% material properties
% anonymous functions f(T) for materials used
% k     [W m-1 K-1]
% c_p   [J kg-1 K-1]
% rho   [kg m-3]
Fk_Cu       = @(T) 385;
Fc_p_Cu     = @(T) 399.5814286 - 0.0585714*T;
Frho_Cu     = @(T) 8940;
Fk_Inco     = @(T) 9.5164989 + 0.0216787*T + -0.0000039*T^2;
Fc_p_Inco   = @(T) 363.8195515 + 0.1233661*T + 0.0000527*T^2;
Frho_Inco   = @(T) 8190;
Fk_H25      = @(T) 9.9905357 + 0.0205437*T - 0.000003*T^2;
Fc_p_H25    = @(T) 396.5228931 + 0.2075422*T + 0.0000134*T^2;
Frho_H25    = @(T) 9070;

% upstream
% rows: k, c_p, rho
% columns: material 1, material 2
mu = cell(3, 2);
mu{1, 1} = Fk_Cu;   mu{2, 1} = Fc_p_Cu;     mu{3, 1} = Frho_Cu;
mu{1, 2} = Fk_Inco; mu{2, 2} = Fc_p_Inco;   mu{3, 2} = Frho_Inco;

% central
% rows: k, c_p, rho
% columns: material 1, material 2, material 3
mc = cell(3, 3);
mc{1, 1} = Fk_Inco; mc{2, 1} = Fc_p_Inco;   mc{3, 1} = Frho_Inco;
mc{1, 2} = Fk_H25;  mc{2, 2} = Fc_p_H25;     mc{3, 2} = Frho_H25;
mc{1, 3} = Fk_Inco; mc{2, 3} = Fc_p_Inco;   mc{3, 3} = Frho_Inco;

% downstream
% rows: k, c_p, rho
% columns: material 1, material 2
md = cell(3, 2);
md{1, 1} = Fk_Inco; md{2, 1} = Fc_p_Inco;   md{3, 1} = Frho_Inco;
md{1, 2} = Fk_Cu;   md{2, 2} = Fc_p_Cu;     md{3, 2} = Frho_Cu;


%% parameters
dt = .25e-3;    % finite volume time step

r = 5; % future points of data considered
epsilon = 1e-2; % small fraction for newton method
RTOLq = 1e-3;    % tolerance for relative change in q
RTOLdq = 1e-2;  % tolerance for relative change in dq
RTOLerror = 1e-6;   % tolerance for relative change in error
maxIter = 20;   % max iterations
qInitial = 10000;   % first guess at q


%% import
dat = readtable(filename);
totalSteps = length(dat.time);
qu = zeros(totalSteps-r-1, 2);
qd = qu;
qu(1, :) = qInitial;
qd(1, :) = qInitial;


%% time iteration upstream
% initial temp distribution
initialT = interp1([1; Nu], [dat.T_Cu2(1); dat.T_Inco1(1)], 1:Nu)';

for m = 1:totalSteps-r-1
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
    T = zeros(Nu, steps(end)+1);
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
            T(:, n+1) = temp1I_PC_QBC(qL, qR, Tavg, T(:, n), dt, gu, ...
                mu, Au);
            TdqL(:, n+1) = temp1I_PC_QBC((1+epsilon)*qL, qR, Tavg, ...
                TdqL(:, n), dt, gu, mu, Au);
            TdqR(:, n+1) = temp1I_PC_QBC(qL, (1+epsilon)*qR, Tavg, ...
                TdqR(:, n), dt, gu, mu, Au);
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
        initialT = temp1I_PC_QBC(qL, qR, Tavg, initialT, dt, gu, mu, Au);
    end

    fprintf('%.3f%% complete.\n', 100*m/(totalSteps-r-1))
end

%% time iteration downstream
% initial temp distribution
initialT = interp1([1; Nd], [dat.T_Inco2(1); dat.T_Cu3(1)], 1:Nd)';

for m = 1:totalSteps-r-1
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
    T = zeros(Nd, steps(end)+1);
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
            T(:, n+1) = temp1I_PC_QBC(qL, qR, Tavg, T(:, n), dt, gd, ...
                md, Ad);
            TdqL(:, n+1) = temp1I_PC_QBC((1+epsilon)*qL, qR, Tavg, ...
                TdqL(:, n), dt, gd, md, Ad);
            TdqR(:, n+1) = temp1I_PC_QBC(qL, (1+epsilon)*qR, Tavg, ...
                TdqR(:, n), dt, gd, md, Ad);
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
        initialT = temp1I_PC_QBC(qL, qR, Tavg, initialT, dt, gd, md, Ad);
    end

    fprintf('%.3f%% complete.\n', 100*m/(totalSteps-r-1))
end

qu = qu(:, 2);
qd = qd(:, 1);
newfilename = filename+".mat";
save(newfilename, "qu", "qd")

