function [datOut] = importInterpData(filename, dt)
% import data from csv & linearly interpolate in time
    % read in data
    datIn = readtable(filename);
    startTime = datIn.time(1);
    endTime = datIn.time(end);
    % vector of times 0-end spaced by dt
    datOut.time = startTime:dt:endTime;
    % interpolate
    datOut.T_Cu1 = interp1(datIn.time, datIn.T_Cu1, datOut.time);
    datOut.T_Cu2 = interp1(datIn.time, datIn.T_Cu2, datOut.time);
    datOut.T_Inco1 = interp1(datIn.time, datIn.T_Inco1, datOut.time);
    datOut.T_Inco2 = interp1(datIn.time, datIn.T_Inco2, datOut.time);
    datOut.T_Cu3 = interp1(datIn.time, datIn.T_Cu3, datOut.time);
    datOut.T_Cu4 = interp1(datIn.time, datIn.T_Cu4, datOut.time);
end

