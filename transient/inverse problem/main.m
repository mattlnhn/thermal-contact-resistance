% %% PROGRESS BAR
% bar = waitbar(0, 'Initialising');

%% MATERIAL PROPERTIES
% anonymous functions so that properties can vary with T
% material 0: Copper
mat.Fk0 = @(T) 385; % []
mat.Fc_p0 = @(T) 399.5814286 - 0.0585714*T; % []
mat.Frho0 = @(T) 8940; % []

% material 1: Inconel 718
mat.Fk1 = @(T) 9.5164989 + 0.0216787*T + -0.0000039*T^2; % []
mat.Fc_p1 = @(T) 363.8195515 + 0.1233661*T + 0.0000527*T^2; % []
mat.Frho1 = @(T) 8190; % []

% material 2: Haynes 25
mat.Fk2 = @(T) 9.9905357 + 0.0205437*T - 0.000003*T^2; % []
mat.Fc_p2 = @(T) 396.5228931 + 0.2075422*T + 0.0000134*T^2; % []
mat.Frho2 = @(T) 9070; % []

% material 3: Inconel 718
mat.Fk3 = @(T) 9.5164989 + 0.0216787*T + -0.0000039*T^2; % []
mat.Fc_p3 = @(T) 363.8195515 + 0.1233661*T + 0.0000527*T^2; % []
mat.Frho3 = @(T) 8190; % []

H = .7875e3;

%% GEOMETRY
% upstream section
geom.L1_u = 1.5e-3; % [m]
geom.L2_u = 6.5e-3; % [m]
geom.dx_u = .25e-3; % [m]

% downstream section
geom.L1_d = 1.5e-3; % [m]
geom.L2_d = 5.0e-3; % [m]
geom.L3_d = 1.5e-3; % [m]
geom.dx_d = 1.0e-4; % [m]

% node positions, referenced to x=0 at TC0
N_u = (geom.L1_u + geom.L2_u)/geom.dx_u;
N_d = (geom.L1_d + geom.L2_d + geom.L3_d)/geom.dx_d;
x_u = ((1:N_u) - .5)*geom.dx_u;
x_d = 8e-3 + ((1:N_d) - .5)*geom.dx_d;
% format for regression
x_u = [ones(N_u, 1) x_u'];
x_d = [ones(N_d, 1) x_d'];
% proportion of section to use for gradient
p_u = ceil(N_u*3/16):N_u;
p_d = 1:floor(N_d*3/16);

%% IMPORT DATA
% read from file
dat = readtable('tcc_data.csv');
[TC0, w0] = smoothdata(dat.T_Cu2, 1, "gaussian");
[TC1, w1] = smoothdata(dat.T_Inco1, 1, "gaussian");
[TC2, w2] = smoothdata(dat.T_Inco2, 1, "gaussian");
time = dat.time;
clear dat

% select range
startTime = 10000;
endTime = startTime + 120;
TC0 = TC0(time>startTime & time<endTime);
TC1 = TC1(time>startTime & time<endTime);
TC2 = TC2(time>startTime & time<endTime);
time = time(time>startTime & time<endTime);

%% TIME ITERATION
dt = .25e-3; % [s]
timeSteps = round((time(end)-time(1))/dt);
d = 1; % data counter

Ti_u = .5*(TC0(1)+TC1(1))*ones(N_u, 1); % initialise with uniform avg temp
Ti_d = .5*(TC1(1)+TC2(1))*ones(N_d, 1);

Tstore_u = zeros(N_u, timeSteps); % initialise storage matrices
Tstore_d = zeros(N_d, timeSteps);
gradients = zeros(2, timeSteps);

for i = 1:timeSteps
    % update boundary conditions
    currentTime = round((i-1)*dt+time(1), 5);
    if currentTime == round(time(d+1), 5)
        d = d+1;
        progress = i/timeSteps;
%         waitmsg = sprintf('%.2f%% complete.', 100*progress);
%         waitbar(progress, bar, waitmsg);
        fprintf('%.2f%% complete\n', 100*progress)
    end

    Ti_u = oneInterfaceNoCR(TC0(d), TC1(d), Ti_u, dt, geom, mat);
    Ti_d = twoInterfaceCR(TC1(d), TC2(d), Ti_d, dt, H, geom, mat);

    beta_u = x_u(p_u, :)\Ti_u(p_u);
    beta_d = x_d(p_d, :)\Ti_d(p_d);
    gradients(:, i) = [beta_u(2) beta_d(2)];

    Tstore_u(:, i) = Ti_u;
    Tstore_d(:, i) = Ti_d;
end

figure(); hold on; plot(gradients(1, :), 'r'); plot(gradients(2, :), 'b')
legend({'upstream temp gradient', 'downstream temp gradient'});
ylim([-480 -440]);
 

