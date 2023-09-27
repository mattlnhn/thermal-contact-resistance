% filename
filename = "230927-tophatsep5-noise1e-04.dat";

% geometry
geometry = cell(6, 1);
geometry{1} = [41.5e-3 8e-3 5e-3 8e-3 14.5e-3];
geometry{2} = sum(geometry{1}, "all");
geometry{3} = .25e-3; % dx
geometry{4} = round(geometry{1}./geometry{3});
geometry{5} = sum(geometry{4});
% positions = [40e-3 44e-3 48e-3 51e-3 53e-3 56e-3 60e-3 64e-3];
positions = [40e-3 48e-3 52e-3 56e-3 64e-3];
geometry{6} = round(.5+(positions)/geometry{3});

% materials
materials = cat(3, materiallookup("Copper"), ...
    materiallookup("Inconel 718"), ...
    materiallookup("Haynes 25"), ...
    materiallookup("Inconel 718"), ...
    materiallookup("Copper"));

% parameters
parameters.r = 2;
parameters.dt = .25e-3;
parameters.epsilon = 1e-3;
parameters.hInitial = [1000 3000];
parameters.maxiter = 10;
parameters.RTOLh = 1e-3;
parameters.RTOLerror = 1e-6;

inverse_sepextra(filename, geometry, materials, parameters)
