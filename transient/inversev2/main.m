%% Sequential function specification method for IHCP with unknown H
% Matt Lenahan, August 2023
% Prerequisites
% 1.    place csv containing temp data in same directory as this file. see
%       example file for formatting and column names.
% 2.    precalculate heat flux in upstream/downstream section using
%       upstreamdownstream.m.
% Method
% The sequential function specification method assumes that at each time
% step, the HTC r steps into the future does not vary appreciably. At each
% time step m: the routine starts with the value of H from the prev. time
% step; computes the temp. distribution with htc = H and with htc = H+eH, 
% where e is a small number, for t = m to (m+r-1); computes the sensitivity
% of temp. to change in H; computes step size dH based on the difference
% between measured and calculated temp. and the sensitivity coefficients;
% and updates H to H+dH and checks for convergence, before moving on to the
% next time step.
% The heat fluxes at the left & right boundaries are prescribed by heat
% flux in the neighbouring sections, calculated from the measured temp. and
% known properties, and the algorithm attempts to match measured & 
% calculated temperature at the boundaries by finding H.

%% prep
clear; close all; clc;
filename = 'data-16-08-2023-12-00.csv';

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

%% precaulculated upstream/downstream data
load upstreamdownstream2.mat

%% central section
dt = 1.25e-3;
tfactor = dt/dt_ud;
dat = importInterpData(filename, dt, 'makima');
timeSteps = length(dat.time);

% initial temp distribution
% linear interpolation between initial end temps
Tc = interp1([1; gc.N], [dat.T_Inco1(1); dat.T_Inco2(1)], 1:gc.N)';

% no. of future time steps
r = 150;
% small number
epsilon = 1e-6; % was 1e-2

% unknown heat transfer coefficient
initialH = 1e3;
Hstore = zeros(timeSteps, 1);

% for calculations
Tcurrent = zeros(gc.N, r);
TdHcurrent = zeros(gc.N, r);
T = zeros(2, r);
TdH = zeros(2, r);
X = zeros(2, r);

% convergence criteria
TOLdH = 1.0e-8;
TOLerror = 1.0e-8; 

% resize heat flux vectors
qu = qu(1:tfactor:end);
qd = qd(1:tfactor:end);

for m = 1:timeSteps % time step m
    if m ~= 1
        H = Hstore(m-1);
    else
        H = initialH;
    end

    % initialise
    converged = 0;
    iter = 1;
    error = 1e6;
    dH = 1e6;
    
    % measured data
    Y = [dat.T_Inco1(m:m+r-1); dat.T_Inco2(m:m+r-1)]; % [2 x r]
    Tavg = .5*(Y(1, :) + Y(2, :));

    while converged == 0
        % previous values for convergence criteria
        prevdH = dH;
        preverror = error;
        
        eH = epsilon*H;

        % starting from prev timestep
        Tcurrent(:, 1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(1), Tc, dt, H, gc, mc);
        TdHcurrent(:, 1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(1), Tc, dt, H+eH, gc, mc);
        for j = 2:r
            % update to (m+j-1)th time step
            Tcurrent(:, j) = temp2I_CR_QBC(qu(m+j-1), qd(m+j-1), Tavg(j), Tcurrent(:, j-1), dt, H, gc, mc);
            TdHcurrent(:, j) = temp2I_CR_QBC(qu(m+j-1), qd(m+j-1), Tavg(j), TdHcurrent(:, j-1), dt, H+eH, gc, mc);
        end
        T = [Tcurrent(1, :); Tcurrent(end, :)];
        TdH = [TdHcurrent(1, :); TdHcurrent(end, :)];
        X = (TdH - T)./eH;

        % newton method update step
        dH = sum((Y-T).*X, 'all')/sum(X.^2, 'all');
        H = H + dH;

        % convergence checks
        error = sum(.5*(Y-T).^2);
        if abs((dH-prevdH)/prevdH) < TOLdH
            converged = 1;
            fprintf('TOLdH. ')
        end
        if abs((error-preverror)/preverror) < TOLerror
            converged = 1;
            fprintf('TOLerror. ')
        end
        if iter > 100
            converged = 1;
            fprintf('Max iter. ')
        end

        iter = iter + 1;

    end

    % to deal with startup error
    if H < 0
        Hstore(m) = 1;
        fprintf('H < 0 caught. ')
    else
        Hstore(m) = H;
    end

    Tc = Tcurrent(:, 1);
    fprintf('%.3f%% complete.\n', 100*m/timeSteps)
end

fileTime = string(datetime('now', Format='yyyyMMdd_HHmmss'));
writematrix(Hstore, 'H_'+fileTime+'.csv');

