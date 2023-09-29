TL = 150;
TR = 25;
totalTime = 180;
warmupTime = 90;
dt = .25e-3; 
DT = 1; % measurement times
ft = DT/dt;
dx = .25e-3;
totalSteps = totalTime/dt;
warmupSteps = warmupTime/dt;
totalSTEPS = totalTime/DT;
% h = interp1([1, .5*(1+totalSteps), totalSteps], [1000 5000 1000], 1:totalSteps, "linear");
h = interp1([0, totalSteps/2, totalSteps], [1000 5000 1000], 1:totalSteps, "nearest");
h0 = 1000;
% t = 1:totalSteps;
% h = 3000 + 2000*(t/totalSteps).*sin(t*4*pi/totalSteps);
% h0 = 3000;

geometry = cell(6, 1);
geometry{1} = [41.5e-3 8e-3 5e-3 8e-3 14.5e-3];
geometry{2} = sum(geometry{1}, "all");
geometry{3} = dx;
geometry{4} = round(geometry{1}./geometry{3});
geometry{5} = sum(geometry{4});
positions = [40e-3 48e-3 56e-3 64e-3];
geometry{6} = round(.5+(positions)/geometry{3});

materials = cat(3, materiallookup("Copper"), ...
                materiallookup("Inconel 718"), ...
                materiallookup("Haynes 25"), ...
                materiallookup("Inconel 718"), ...
                materiallookup("Copper"));

N = geometry{5};
A = diag(ones(N-1, 1), -1) + diag(-2*ones(N, 1), 0) + ...
    diag(ones(N-1, 1), 1);
A(1, 1) = -3; A(end, end) = -3;

Tsave = zeros(N, totalSTEPS);
T = interp1([1 N]', [TR TR]', (1:N)');

%% warmup
for w = 1:warmupSteps
    Tavg = .5*([TL; T(geometry{6})]+[T(geometry{6}); TR]);
    T = direct([TL TR], Tavg, T, h0, dt, geometry, materials, A);
end

%% generate
for m = 1:totalSteps
    %Tavg = ones(5, 1)*mean(T, "all");
    Tavg = .5*([TL; T(geometry{6})]+[T(geometry{6}); TR]);
    T = direct([TL TR], Tavg, T, h(m), dt, geometry, materials, A);
    if mod(m, ft) == 0 % every big time step elapsed
        Tsave(:, m/ft) = T;
    end
end

Tout = [TL*ones(totalSTEPS, 1), Tsave(geometry{6}, :)', TR*ones(totalSTEPS, 1)];
[rows, cols] = size(Tout);
noise = 1e-3;
accuracy = 1.5;

noisyTout = Tout + noise*Tout.*(rand(rows, cols) - .5) + ...
    accuracy*ones(rows, 1).*(rand(1, cols) - .5); % rnd noise + constant offset
Touttable = array2table([(DT:DT:totalTime)' noisyTout]);
Touttable.Properties.VariableNames(1:7) = ["time", "T_Cu1", "T_Cu2", "T_Inco1", "T_Inco2", "T_Cu3", "T_Cu4"];
filename = sprintf("230926-tophat-noise%.0e-pm%.1f", noise, accuracy) + ".dat";
writetable(Touttable, filename)


