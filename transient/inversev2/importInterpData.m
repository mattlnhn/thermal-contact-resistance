function [datOut] = importInterpData(filename, dt, method, smoothing, window)
% import data from csv & linearly interpolate in time
    % read in data
    datIn = readtable(filename);
    startTime = datIn.time(1);
    endTime = datIn.time(end);
    % smooth temps before interp
    if exist('smoothing', 'var') && exist('window', 'var')
        datIn.T_Cu1 = smoothdata(datIn.T_Cu1, smoothing, window);
        datIn.T_Cu2 = smoothdata(datIn.T_Cu2, smoothing, window);
        datIn.T_Inco1 = smoothdata(datIn.T_Inco1, smoothing, window);
        datIn.T_Inco2 = smoothdata(datIn.T_Inco2, smoothing, window);
        datIn.T_Cu3 = smoothdata(datIn.T_Cu3, smoothing, window);
        datIn.T_Cu4 = smoothdata(datIn.T_Cu4, smoothing, window);
    end
    % vector of times 0-end spaced by dt
    datOut.time = startTime:dt:endTime;
    % interpolate
    datOut.load = interp1(datIn.time, datIn.load, datOut.time, method);
    datOut.T_Cu1 = interp1(datIn.time, datIn.T_Cu1, datOut.time, method);
    datOut.T_Cu2 = interp1(datIn.time, datIn.T_Cu2, datOut.time, method);
    datOut.T_Inco1 = interp1(datIn.time, datIn.T_Inco1, datOut.time, method);
    datOut.T_Inco2 = interp1(datIn.time, datIn.T_Inco2, datOut.time, method);
    datOut.T_Cu3 = interp1(datIn.time, datIn.T_Cu3, datOut.time, method);
    datOut.T_Cu4 = interp1(datIn.time, datIn.T_Cu4, datOut.time, method);
end

