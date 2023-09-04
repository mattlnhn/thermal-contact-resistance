%% prep
clear;

%% plot
figure()

subplot(1, 2, 1)
title(["5N loading/unloading cycles, 150 deg C", "New algorithm, r = 5"])
xlabel('Pressure [kPa]')
ylabel('Contact conductance [m^2 K W^{-1}]')
hold on
grid minor
ylim([0 1e-3])

subplot(1, 2, 2)
xlabel('Time [s]')
hold on
grid minor
yyaxis left
ylabel('HTC [W m^{-2} K^{-1}]')
ylim([1000 4000])
yyaxis right
ylabel('Pressure [kPa]')
ylim([0 25])

cmap = hsv(8);


% dat = readtable("20230830_1_5Ncycle1.dat");
% pressure = -dat.load*4.8217; % kPa
% time = dat.time;
% h = readmatrix("h2_20230830_1_5Ncycle1.dat");
% r = 5;
% n = 1;
% 
% subplot(1, 2, 1)
% plot(pressure(1:end-r), h(:, 1).^-1, '-', Color=cmap(n, :))
% plot(pressure(1:end-r), h(:, 2).^-1, ':', Color=cmap(n, :))
% subplot(1, 2, 2)
% yyaxis left
% plot(time(1:end-r), h(:, 1), '-', Color=cmap(n, :))
% plot(time(1:end-r), h(:, 2), ':', Color=cmap(n, :))
% yyaxis right
% plot(time, pressure, ':', Color=cmap(n, :))
% 
% 
% dat = readtable("20230830_1_5Ncycle2.dat");
% pressure = -dat.load*4.8217; % kPa
% time = dat.time;
% h = readmatrix("h2_20230830_1_5Ncycle2.dat");
% r = 5;
% n = 2;
% 
% subplot(1, 2, 1)
% plot(pressure(1:end-r), h(:, 1).^-1, '-', Color=cmap(n, :))
% plot(pressure(1:end-r), h(:, 2).^-1, ':', Color=cmap(n, :))
% subplot(1, 2, 2)
% yyaxis left
% plot(time(1:end-r), h(:, 1), '-', Color=cmap(n, :))
% plot(time(1:end-r), h(:, 2), ':', Color=cmap(n, :))
% yyaxis right
% plot(time, pressure, ':', Color=cmap(n, :))
% 
% 
% dat = readtable("20230830_1_5Ncycle3.dat");
% pressure = -dat.load*4.8217; % kPa
% time = dat.time;
% h = readmatrix("h2_20230830_1_5Ncycle3.dat");
% r = 5;
% n = 3;
% 
% subplot(1, 2, 1)
% plot(pressure(1:end-r), h(:, 1).^-1, '-', Color=cmap(n, :))
% plot(pressure(1:end-r), h(:, 2).^-1, ':', Color=cmap(n, :))
% subplot(1, 2, 2)
% yyaxis left
% plot(time(1:end-r), h(:, 1), '-', Color=cmap(n, :))
% plot(time(1:end-r), h(:, 2), ':', Color=cmap(n, :))
% yyaxis right
% plot(time, pressure, ':', Color=cmap(n, :))
% 
% 
% dat = readtable("20230830_1_5Ncycle4.dat");
% pressure = -dat.load*4.8217; % kPa
% time = dat.time;
% h = readmatrix("h2_20230830_1_5Ncycle4.dat");
% r = 5;
% n = 4;
% 
% subplot(1, 2, 1)
% plot(pressure(1:end-r), h(:, 1).^-1, '-', Color=cmap(n, :))
% plot(pressure(1:end-r), h(:, 2).^-1, ':', Color=cmap(n, :))
% subplot(1, 2, 2)
% yyaxis left
% plot(time(1:end-r), h(:, 1), '-', Color=cmap(n, :))
% plot(time(1:end-r), h(:, 2), ':', Color=cmap(n, :))
% yyaxis right
% plot(time, pressure, ':', Color=cmap(n, :))
% 
% 
% dat = readtable("20230830_1_5Ncycle6.dat");
% pressure = -dat.load*4.8217; % kPa
% time = dat.time;
% h = readmatrix("h2_20230830_1_5Ncycle6.dat");
% r = 5;
% n = 6;
% 
% subplot(1, 2, 1)
% plot(pressure(1:end-r), h(:, 1).^-1, '-', Color=cmap(n, :))
% plot(pressure(1:end-r), h(:, 2).^-1, ':', Color=cmap(n, :))
% subplot(1, 2, 2)
% yyaxis left
% plot(time(1:end-r), h(:, 1), '-', Color=cmap(n, :))
% plot(time(1:end-r), h(:, 2), ':', Color=cmap(n, :))
% yyaxis right
% plot(time, pressure, ':', Color=cmap(n, :))
% 
% 
% dat = readtable("20230830_1_5Ncycle7.dat");
% pressure = -dat.load*4.8217; % kPa
% time = dat.time;
% h = readmatrix("h2_20230830_1_5Ncycle7.dat");
% r = 5;
% n = 7;
% 
% subplot(1, 2, 1)
% plot(pressure(1:end-r), h(:, 1).^-1, '-', Color=cmap(n, :))
% plot(pressure(1:end-r), h(:, 2).^-1, ':', Color=cmap(n, :))
% subplot(1, 2, 2)
% yyaxis left
% plot(time(1:end-r), h(:, 1), '-', Color=cmap(n, :))
% plot(time(1:end-r), h(:, 2), ':', Color=cmap(n, :))
% yyaxis right
% plot(time, pressure, ':', Color=cmap(n, :))


dat = readtable("20230830_1_5Ncycle8.dat");
pressure = -dat.load*4.8217; % kPa
time = dat.time;
h = readmatrix("h2_20230830_1_5Ncycle8.dat");
r = 10;
n = 8;

subplot(1, 2, 1)
plot(pressure(1:end-r), h(:, 1).^-1, '-', Color=cmap(n, :))
plot(pressure(1:end-r), h(:, 2).^-1, ':', Color=cmap(n, :))
subplot(1, 2, 2)
yyaxis left
plot(time(1:end-r), h(:, 1), '-', Color=cmap(n, :))
plot(time(1:end-r), h(:, 2), ':', Color=cmap(n, :))
yyaxis right
plot(time, pressure, ':', Color=cmap(n, :))


