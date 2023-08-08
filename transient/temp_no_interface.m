% infer transient temp distribution between 2 TCs

clear; close all;

% import data, assign to vars & clear duplicates
dat = readtable('tcc_data.csv');
time = dat.time;    % time
TC1 = dat.T_Cu1;    % temp
TC2 = dat.T_Cu2;
TC3 = dat.T_Inco1;
TC4 = dat.T_Inco2;
TC5 = dat.T_Cu3;
TC6 = dat.T_Cu4;
roomTemp = mean(dat.T_room);
clear dat

% geometry
dx = .5e-3;
L = .04; % distance between TCs, m

% time
dt = 1e-3;
timeSteps = time(end)/dt;

% copper section 1, TC1-TC2
% functions for varying physical properties w/ temp
Fk = @(T) 385;
Fc_p = @(T) 399.5814286 - 0.0585714*T;
Frho = @(T) 8940;

% nodes
N = L/dx;

% constructing matrices
e = ones(N, 1);
A = spdiags([e -2*e e], [-1 0 1], N, N);
A = full(A);
A(1, 1) = -3; A(end, end) = -3;

b = zeros(N, 1);

Ti = roomTemp*ones(N, 1); % assume room temp throughout at start
Tstore = zeros(N, timeSteps);

d = 1; % counts through data points

for i = 1:timeSteps
    % update boundary conditions when next data point reached
    currentTime = round((i-1)*dt, 5);
    if currentTime == round(time(d+1), 5)
        d = d+1;
        fprintf('%.2f%% complete.\n', 100*i/timeSteps)
    end
    b(1) = 2*TC1(d);
    b(end) = 2*TC2(d);

    % physical properties based on average temp
    Tavg = .5*(TC1(d) + TC2(d));
    k = Fk(Tavg);
    c_p = Fc_p(Tavg);
    rho = Frho(Tavg);
    alpha = k*rho^-1*c_p^-1;
    tau = dt*alpha*dx^-2;

    % break if stability criterion not met
    if tau > .5
        fprintf('Unstable. tau = %d', tau)
        break
    end

    % T at next time step
    Ti = tau*(A*Ti + b) + Ti;

    % write current state to memory
    Tstore(:, i) = Ti;
end

