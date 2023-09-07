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
filename = "Data Record -06-09-2023- 14-49_tx.dat";


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
gc = cell(3, 1);
gc{1} = 1.0e-4;
gc{2} = [1.5e-3, 5.0e-3, 1.5e-3];
gc{3} = [gc{2}(1)/gc{1}, gc{2}(2)/gc{1}, gc{2}(3)/gc{1}];
Nc = round(sum(gc{3}, "all"));
TCH25 = floor(Nc*2.5e-3/8e-3); % node # for H25 thermocouple
% downstream inco/copper
gd = cell(3, 1);
gd{1} = .25e-3;
gd{2} = [6.5e-3, 1.5e-3];
gd{3} = [gd{2}(1)/gd{1}, gd{2}(2)/gd{1}];
Nd = round(sum(gd{3}, "all"));


%% precalculate A matrices
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
% coefficients, f(T) = a + bT + cT^2
k_Cu = [385 0 0];
c_p_Cu = [399.5814286 -.0585714 0];
rho_Cu = [8940 0 0];

k_Inco = [9.5164989 .0216787 -.0000039];
c_p_Inco = [363.8195515 .1233661 .0000527];
rho_Inco = [8190 0 0];

k_H25 = [9.9905357 .0205437 -.000003];
c_p_H25 = [396.5228931 .2075422 .0000134];
rho_H25 = [9070 0 0];

% upstream
% rows: k, c_p, rho
% columns: a, b, c where f(T) = a + bT + cT^2
% layers: material
mu = zeros(3, 3, 2);
mu(1, :, 1) = k_Cu;     mu(2, :, 1) = c_p_Cu;   mu(3, :, 1) = rho_Cu;
mu(1, :, 2) = k_Inco;   mu(2, :, 2) = c_p_Inco; mu(3, :, 2) = rho_Inco;

% central
mc = zeros(3, 3, 3);
mc(1, :, 1) = k_Inco;   mc(2, :, 1) = c_p_Inco; mc(3, :, 1) = rho_Inco;
mc(1, :, 2) = k_H25;    mc(2, :, 2) = c_p_H25;  mc(3, :, 2) = rho_H25;
mc(1, :, 3) = k_Inco;   mc(2, :, 3) = c_p_Inco; mc(3, :, 3) = rho_Inco;

% downstream
md = zeros(3, 3, 2);
md(1, :, 1) = k_Inco;   md(2, :, 1) = c_p_Inco; md(3, :, 1) = rho_Inco;
md(1, :, 2) = k_Cu;     md(2, :, 2) = c_p_Cu;   md(3, :, 2) = rho_Cu;


%% precaulculated upstream/downstream data
load(filename+".mat")


%% parameters
dt = .25e-3;    % finite volume time step

r = 20; % future points of data considered
epsilon = 1e-6; % small fraction for newton method
RTOLh = 1e-3;    % tolerance for relative change in h
RTOLerror = 1e-6;   % tolerance for relative change in error
maxIter = 20;   % max iterations    
hInitial = 1000;   % first guess at h


%% import
dat = readtable(filename);
totalSteps = length(dat.time);
hStore = zeros(totalSteps-r-1, 2);
hStore(1, :) = hInitial;
figure(); hold on;
plot(dat.time, dat.T_oven, 'k--')
plot(dat.time, dat.T_Cu1, 'r-')
plot(dat.time, dat.T_Cu2, 'r--')
plot(dat.time, dat.T_Inco1, 'b-')
plot(dat.time, dat.T_H25, 'g-')
plot(dat.time, dat.T_Inco2, 'b--')
plot(dat.time, dat.T_Cu3, 'r:')


%% time iteration
% initial temp distribution
initialT = interp1([1; TCH25; Nc], [dat.T_Inco1(1); dat.T_H25(1); dat.T_Inco2(1)], (1:Nc)');

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
    % no. of steps of dt from initial time
    steps = round((time - time(1))/dt);

    % temp matrices incl. initial time
    T = zeros(Nc, steps(end)+1);
    T(:, 1) = initialT;
    Tdhu = T;
    Tdhd = T;

    % start with previous values of h
    hu = hStore(m, 1);
    hd = hStore(m, 2);

    % measured data r steps into future
    Y = [dat.T_Inco1(m:m+r)';
        dat.T_H25(m:m+r)';
        dat.T_Inco2(m:m+r)'];
    Tavg = mean(Y, 1);

    while converged == 0
        prevdhu = dhu;
        prevdhd = dhd;
        preverroru = erroru;
        preverrord = errord;

        for n = 1:steps(end)
            % index for avg temp
            ai = sum(steps<n);
            % temp updates
            T(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                T(:, n), dt, [hu hd], gc, mc, Ac);
            Tdhu(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                Tdhu(:, n), dt, [(1+epsilon)*hu hd], gc, mc, Ac);
            Tdhd(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                Tdhd(:, n), dt, [hu (1+epsilon)*hd], gc, mc, Ac);
        end

        % T at sensor locations & measurement times, excl. initial
        Ts = [T(1, steps+1); T(TCH25, steps+1); T(end, steps+1)];
        Tdhus = [Tdhu(1, steps+1); Tdhu(TCH25, steps+1); Tdhu(end, steps+1)];
        Tdhds = [Tdhd(1, steps+1); Tdhd(TCH25, steps+1); Tdhd(end, steps+1)];

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
        initialT = temp2I_CR_QBC(qu(m), qd(m), Tavg(1), ...
                initialT, dt, [hu hd], gc, mc, Ac);
    end

    fprintf('%.3f%% complete.\n', 100*m/(totalSteps-r-1))
end

writematrix(hStore, "h2_"+filename);

