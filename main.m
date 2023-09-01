%% Function specification method for IHCP with unknown H
% Matt Lenahan, August 2023
% Prerequisites
% 1.    place data file containing temp data in same directory as this
%       file. variable (column) names are:
%       pos     LVDT position
%       load    load
%       time    time (seconds of the day)
%       T_oven, T_Cu1, T_Cu2, T_Inco1, T_H25, T_Inco2, T_Cu3, T_Cu4
%               thermocouples, in order from oven side to load cell side
% 2.    precalculate heat flux in upstream/downstream section using
%       upstreamdownstream.m.

%% prep
clear;
filename = "20230831_cycles_fixed_2of3.txt";

%% section geometry
% upstream copper/inco
gu.dx = .25e-3;
gu.L1 = 1.5e-3;
gu.L2 = 6.5e-3;
gu.N = (gu.L1+gu.L2)/gu.dx;
% central inco/h25/inco
gc.dx = 1.0e-4;
gc.L1 = 1.5e-3;
gc.L2 = 5.0e-3;
gc.L3 = 1.5e-3;
gc.N = (gc.L1+gc.L2+gc.L3)/gc.dx;
gc.tch25 = floor(gc.N*1.5e-3/5e-3);
% downstream inco/copper
gd.dx = .25e-3;
gd.L1 = 6.5e-3;
gd.L2 = 1.5e-3;
gd.N = (gd.L1+gd.L2)/gd.dx;

%% material properties
% anonymous functions f(T) for materials used
% k     [W m-1 K-1]
% c_p   [J kg-1 K-1]
% rho   [kg m-3]
Fk_Cu = @(T) 385;
Fc_p_Cu = @(T) 399.5814286 - 0.0585714*T;
Frho_Cu = @(T) 8940;
Fk_Inco = @(T) 9.5164989 + 0.0216787*T + -0.0000039*T^2;
Fc_p_Inco = @(T) 363.8195515 + 0.1233661*T + 0.0000527*T^2;
Frho_Inco = @(T) 8190;
Fk_H25 = @(T) 9.9905357 + 0.0205437*T - 0.000003*T^2;
Fc_p_H25 = @(T) 396.5228931 + 0.2075422*T + 0.0000134*T^2;
Frho_H25 = @(T) 9070;

% upstream
mu.Fk1 = Fk_Cu; mu.Fc_p1 = Fc_p_Cu; mu.Frho1 = Frho_Cu;
mu.Fk2 = Fk_Inco; mu.Fc_p2 = Fc_p_Inco; mu.Frho2 = Frho_Inco;

% central
mc.Fk1 = Fk_Inco; mc.Fc_p1 = Fc_p_Inco; mc.Frho1 = Frho_Inco;
mc.Fk2 = Fk_H25; mc.Fc_p2 = Fc_p_H25; mc.Frho2 = Frho_H25;
mc.Fk3 = Fk_Inco; mc.Fc_p3 = Fc_p_Inco; mc.Frho3 = Frho_Inco;

% downstream
md.Fk1 = Fk_Inco; md.Fc_p1 = Fc_p_Inco; md.Frho1 = Frho_Inco;
md.Fk2 = Fk_Cu; md.Fc_p2 = Fc_p_Cu; md.Frho2 = Frho_Cu;

%% precaulculated upstream/downstream data
load(filename+".mat")

%% import
dat = readtable(filename);
totalSteps = length(dat.time);
hStore = zeros(totalSteps+1, 2);

%% parameters
dt = 1.25e-3;    % finite volume time step

r = 5; % future points of data considered
epsilon = 1e-2; % small fraction for newton method
RTOLh = 1e-3;    % tolerance for relative change in q
RTOLdh = 1e-2;  % tolerance for relative change in dq
RTOLerror = 1e-6;   % tolerance for relative change in error
maxIter = 20;   % max iterations    

hInitial = 1000;   % first guess at h
hStore(1, :) = hInitial;

%% time iteration
% initial temp distribution
initialT = interp1([1; gc.tch25; gc.N], [dat.T_Inco1(1); dat.T_H25(1); dat.T_Inco2(1)], 1:gc.N)';

for m = 1:totalSteps-r-1
    % reset counters
    converged = 0;
    iter = 1;
    erroru = 1e20;
    errord = erroru;
    dhu = erroru;
    dhd = erroru;

    % time data
    time = dat.time(m:m+r);
    steps = round((time - time(1))/dt);

    % temp matrices incl. initial time
    T = zeros(gc.N, steps(end)+1);
    T(:, 1) = initialT;
    Tdhu = T;
    Tdhd = T;

    hu = hStore(m, 1);
    hd = hStore(m, 2);

    Y = [dat.T_Inco1(m:m+r)';
        dat.T_H25(m:m+r)';
        dat.T_Inco2(m:m+r)'];
    Tavg = mean(Y, 'all');

    while converged == 0
        prevdhu = dhu;
        prevdhd = dhd;
        preverroru = erroru;
        preverrord = errord;

        for n = 1:steps(end)
            T(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg, ...
                T(:, n), dt, [hu hd], gc, mc);
            Tdhu(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg, ...
                Tdhu(:, n), dt, [(1+epsilon)*hu hd], gc, mc);
            Tdhd(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg, ...
                Tdhd(:, n), dt, [hu (1+epsilon)*hd], gc, mc);
        end

        % T at sensor locations & measurement times, excl. initial
        Ts = [T(1, steps+1); T(gc.tch25, steps+1); T(end, steps+1)];
        Tdhus = [Tdhu(1, steps+1); Tdhu(gc.tch25, steps+1); Tdhu(end, steps+1)];
        Tdhds = [Tdhd(1, steps+1); Tdhd(gc.tch25, steps+1); Tdhd(end, steps+1)];

        % sensitivity coefficients
        Xu = (Tdhus - Ts)./(epsilon*hu);
        Xd = (Tdhds - Ts)./(epsilon*hd);

        % newton method increment
        dhu = sum((Y-Ts).*Xu, 'all')/sum(Xu.^2, 'all');
        dhd = sum((Y-Ts).*Xd, 'all')/sum(Xd.^2, 'all');

        hu = hu + dhu;
        hd = hd + dhd;

        erroru = sum(.5*(Y-Tdhus).^2, 'all');
        errord = sum(.5*(Y-Tdhds).^2, 'all');

        % convergence criteria
        % relative change in h
        if abs(dhu/hu) < RTOLh && abs(dhd/hd) < RTOLh
            converged = 1;
            fprintf('RTOLh. ')
        end
        % relative change in error
        if abs((erroru-preverroru)/preverroru) < RTOLerror && ...
                abs((errord-preverrord)/preverrord) < RTOLerror
            converged = 1;
            fprintf('RTOLerror. ')
        end
        % timeout
        if iter > maxIter
            converged = 1;
            fprintf('Max iter. ')
        end

        iter = iter + 1;

    end

    % save q
    hStore(m+1, 1) = hu;
    hStore(m+1, 2) = hd;
    
    % next step initial temp
    for n = 1:steps(2)
        initialT = temp2I_CR_QBC(qu(m), qd(m), Tavg, ...
                initialT, dt, [hu hd], gc, mc);
    end

    fprintf('%.3f%% complete.\n', 100*m/(totalSteps-r-1))
end

writematrix(hStore, "h2_"+filename);

