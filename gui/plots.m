filenames = ["230915-1414.dat", "230915-1441.dat", "230915-1502.dat"];
hfilenames = "h_"+filenames;

% set up figure
figure()
hold on
grid on
grid minor
xlabel("Pressure [Pa]")
ylabel("Contact resistance h^{-1} [m^2 K W^{-1}]")
ylim([0 1e-3])

cmap = hsv(length(filenames));

% plot
r = 2;

for f = 1:length(filenames)
    dat = readtable(filenames(f));
    hdat = readmatrix(hfilenames(f));
    plot(-dat.load(1:end-r)/(16.25e-3^2*pi/4), hdat(2:end).^-1, '-', Color=cmap(f, :))
end