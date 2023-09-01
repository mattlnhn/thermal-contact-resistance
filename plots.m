%% prep
clear;

%% plot
figure()

subplot(1, 2, 1)
title(["5N loading/unloading cycles, 150 deg C", "New algorithm, r = 5"])
xlabel('Pressure [kPa]')
ylabel('HTC [W m^{-2} K^{-1}]')
hold on
grid minor
ylim([1000 4000])

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


dat = readtable("20230831_cycles_fixed_1of3.txt");
pressure = -dat.load*4.8217; % kPa
time = dat.time;
h = readmatrix("h2_20230831_cycles_fixed_1of3.txt");

subplot(1, 2, 1)
plot(pressure, h(2:end, 1), '-', Color=cmap(1, :))
plot(pressure, h(2:end, 2), '--', Color=cmap(1, :))
subplot(1, 2, 2)
yyaxis left
plot(time, h(2:end, 1), '-', Color=cmap(1, :))
plot(time, h(2:end, 2), '--', Color=cmap(1, :))
yyaxis right
plot(time, pressure, ':', Color=cmap(1, :))



dat = readtable("20230831_cycles_fixed_2of3.txt");
pressure = -dat.load*4.8217; % kPa
time = dat.time;
h = readmatrix("h2_20230831_cycles_fixed_2of3.txt");

subplot(1, 2, 1)
plot(pressure, h(2:end, 1), '-', Color=cmap(2, :))
plot(pressure, h(2:end, 2), '--', Color=cmap(2, :))
subplot(1, 2, 2)
yyaxis left
plot(time, h(2:end, 1), '-', Color=cmap(2, :))
plot(time, h(2:end, 2), '--', Color=cmap(2, :))
yyaxis right
plot(time, pressure, ':', Color=cmap(2, :))


dat = readtable("20230831_cycles_fixed_3of3.txt");
pressure = -dat.load*4.8217; % kPa
time = dat.time;
h = readmatrix("h2_20230831_cycles_fixed_3of3.txt");

subplot(1, 2, 1)
plot(pressure, h(2:end, 1), '-', Color=cmap(3, :))
plot(pressure, h(2:end, 2), '--', Color=cmap(3, :))
subplot(1, 2, 2)
yyaxis left
plot(time, h(2:end, 1), '-', Color=cmap(3, :))
plot(time, h(2:end, 2), '--', Color=cmap(3, :))
yyaxis right
plot(time, pressure, ':', Color=cmap(3, :))
