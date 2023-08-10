% infer transient temp distribution between 2 TCs w/ material interface

clear; close all;

% import data, assign to vars & clear duplicates
dat = readtable('tcc_data.csv');
% timestamps
time = dat.time;
% TC temps: TC1 = left, TC2 = right
TC1 = dat.T_Inco1;
TC2 = dat.T_Inco2;
% average room temp
roomTemp = mean(dat.T_room);
clear dat

% truncate data if desired
starttime = 0;
TC1 = TC1(time>starttime);
TC2 = TC2(time>starttime);
time = time(time>starttime);

% geometry
dx = .5e-3; % cell size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L1 = .0015; % span of 1st material, m
L2 = .0050; % span of 2nd material, m
L3 = .0015; % span of 2nd material, m
Ltotal = L1 + L2 + L3; % distance between TCs

% time
dt = 1e-3; % time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeSteps = round((time(end)-time(1))/dt);

% functions for varying physical properties w/ temp
% material 1: Inco
Fk1 = @(T) 9.5164989 + 0.0216787*T + -0.0000039*T^2;
Fc_p1 = @(T) 363.8195515 + 0.1233661*T + 0.0000527*T^2;
Frho1 = @(T) 8190;
% material 2: H25
%Fk2 = @(T) 9.9905357 + 0.0205437*T - 0.000003*T^2;
%Fc_p2 = @(T) 396.5228931 + 0.2075422*T + 0.0000134*T^2;
%Frho2 = @(T) 9070;
% material 2a: copper
Fk2 = @(T) 385;
Fc_p2 = @(T) 399.5814286 - 0.0585714*T;
Frho2 = @(T) 8940;
% material 3: Inco
Fk3 = @(T) 9.5164989 + 0.0216787*T + -0.0000039*T^2;
Fc_p3 = @(T) 363.8195515 + 0.1233661*T + 0.0000527*T^2;
Frho3 = @(T) 8190;

% nodes
N = Ltotal/dx;
N1 = L1/dx;
N2 = L2/dx;
N3 = L3/dx;

% constructing matrices
e = ones(N, 1);
A = spdiags([e -2*e e], [-1 0 1], N, N);
A = full(A);
A(1, 1) = -3; A(end, end) = -3;

b = zeros(N, 1);

%Ti = roomTemp*ones(N, 1); % assume room temp throughout at start
Ti = .5*(TC1(1)+TC2(1))*ones(N, 1); % assume uniform avg temp to start
Tstore = zeros(N, timeSteps);

d = 1; % counts through data points

for i = 1:timeSteps
    % update boundary conditions when next data point reached
    currentTime = round((i-1)*dt+time(1), 5);
    if currentTime == round(time(d+1), 5)
        d = d+1;
        fprintf('%.2f%% complete.\n', 100*i/timeSteps)
    end
    b(1) = 2*TC1(d);
    b(end) = 2*TC2(d);

    % physical properties based on average temp
    Tavg = .5*(TC1(d) + TC2(d));

    k1 = Fk1(Tavg);
    c_p1 = Fc_p1(Tavg);
    rho1 = Frho1(Tavg);
    alpha1 = k1*rho1^-1*c_p1^-1;
    tau1 = dt*alpha1*dx^-2;
    tau1D = tau1*eye(N1); % diag matrix of tau1

    k2 = Fk2(Tavg);
    c_p2 = Fc_p2(Tavg);
    rho2 = Frho2(Tavg);
    alpha2 = k2*rho2^-1*c_p2^-1;
    tau2 = dt*alpha2*dx^-2;
    tau2D = tau2*eye(N2); % diag matrix of tau2

    k3 = Fk3(Tavg);
    c_p3 = Fc_p3(Tavg);
    rho3 = Frho3(Tavg);
    alpha3 = k3*rho3^-1*c_p3^-1;
    tau3 = dt*alpha3*dx^-2;
    tau3D = tau3*eye(N3); % diag matrix of tau3
    
    % break if stability criterion not met
    if tau1 > .5 || tau2 > .5 || tau3 > .5
        fprintf('Unstable. tau1 = %d, tau2 = %d, tau3 = %d\n', tau1, tau2, tau3)
        break
    end
    
    % block diagonal matrix so that diffusivity corresponds to material
    tau = blkdiag(tau1D, tau2D, tau3D);

    % modification of matrix A on boundary
    % see https://engineering.stackexchange.com/questions/5921/modeling-
    % transient-heat-transfer-between-two-1-d-materials
    kappa12 = k1/k2;
    gamma12 = (kappa12-1)/(kappa12+1);
    A(N1, N1) = -(2-gamma12); A(N1+1, N1+1) = -(2+gamma12);
    A(N1, N1+1) = 2/(kappa12+1); A(N1+1, N1) = 2*kappa12/(kappa12+1);

    kappa23 = k2/k3;
    gamma23 = (kappa23-1)/(kappa23+1);
    A(N1+N2, N1+N2) = -(2-gamma23); A(N1+N2+1, N1+N2+1) = -(2+gamma23);
    A(N1+N2, N1+N2+1) = 2/(kappa23+1); A(N1+N2+1, N1+N2) = 2*kappa23/(kappa23+1);

    % T at next time step
    Ti = tau*(A*Ti + b) + Ti;

    % write current state to memory
    Tstore(:, i) = Ti;
end

