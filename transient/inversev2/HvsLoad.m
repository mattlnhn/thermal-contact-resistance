%% prep
clear; close all; clc;
filename_dat = 'data-16-08-2023-12-00.csv';
filename_H = 'H_20230817_104044.csv';

%% import
dt = 1.25e-3;
dat = importInterpData(filename_dat, dt, 'makima');
load = dat.load;
time = dat.time;
clear dat

Hread = readtable('H_20230817_104044.csv');
H = Hread.Var1;

%% truncate
startpt = 5500;
endpt = 200;
time = time(startpt:end-endpt);
load = load(startpt:end-endpt);
H = H(startpt:end-endpt);

%% filter
window = 2000;
method = 'gaussian';
load_filtered = smoothdata(load, method, window);
H_filtered = smoothdata(H, method, window);

%% plots
% H and load vs. time
figure(); hold on;

subplot(2, 2, 1)
plot(time, H, 'r')
grid minor
ylim([850 950])
xlabel('Time, s')
ylabel('H, W m^{-2} K^{-1}')
title('H')

subplot(2, 2, 2)
plot(time, -load, 'b')
grid minor
xlabel('Time, s')
ylabel('Load, N')
title('Load')

subplot(2, 2, 3)
plot(time, H_filtered, 'r')
grid minor
ylim([850 950])
xlabel('Time, s')
ylabel('H, W m^{-2} K^{-1}')
title('H, filtered')

subplot(2, 2, 4)
plot(time, -load_filtered, 'b')
grid minor
xlabel('Time, s')
ylabel('Load, N')
title('Load, filtered')

% H vs load with regression
figure()

% linear regression
X = [ones(length(load), 1) -load'];
y = H;
coef = X\y;
xplot = 0:.1:3;
yplot = coef(1) + coef(2).*xplot;

subplot(1, 2, 1)
hold on
plot(-load, H, 'k:')
plot(xplot, yplot, 'r')
grid minor
ylim([870 970])
xlabel('Load, N')
ylabel('H, W m^{-2} K^{-1}')
title(sprintf('Regression %.3f + %.3f*x', coef(1), coef(2)))

% linear regression
X = [ones(length(load_filtered), 1) -load_filtered'];
y = H_filtered;
coef = X\y;
xplot = 0:.1:3;
yplot = coef(1) + coef(2).*xplot;

subplot(1, 2, 2)
hold on
plot(-load_filtered, H_filtered, 'k:')
plot(xplot, yplot, 'r')
grid minor
ylim([870 970])
xlabel('Load, N')
ylabel('H, W m^{-2} K^{-1}')
title(sprintf('Filtered data, regression %.3f + %.3f*x', coef(1), coef(2)))
