%% TRANSIENT

filenames = ["H:\Internship\MATLAB\git\gui\230928-smooth.dat", ...
    "H:\Internship\MATLAB\git\gui\230928-med.dat", ...
    "H:\Internship\MATLAB\git\gui\230928-rough.dat"];
[path, file, ext] = fileparts(filenames);
hfilenames = path+"\h_"+file+ext;

% set up figure
figure()
subplot(1, 2, 1)
hold on
grid on
grid minor
xlabel("Pressure [Pa]")
ylabel("Contact resistance h^{-1} [m^2 K W^{-1}]")
ylim([0 1e-3])
subplot(1, 2, 2)
hold on
grid on
grid minor
xlabel("Pressure [Pa]")
ylabel("Contact resistance h^{-1} [m^2 K W^{-1}]")
ylim([0 1e-3])

cmap = winter(length(filenames));

% plot
r = 2;

for f = 1:length(filenames)
    dat = readtable(filenames(f));
    hdat = readmatrix(hfilenames(f));
    subplot(1, 2, 1)
    plot(-9.81*dat.load(1:end-r)/(16.25e-3^2*pi/4), hdat(2:end).^-1, '-', Color=cmap(f, :))
end

%% PSEUDO-STEADY

L = length(filenames);

for ff = 1:L
    % data
    dat = readtable(filenames(ff));
    % TC positions referenced to T_Cu1
    TCpos = [0, .04, .048, .056, .064, .077]; % m
    
    % thermal conductivities, W m-1 K-1
    k_Cu = 385; % W m-2 K-1
    func_k_Inco = @(T) 9.5164989 + .0216787*T + -.0000039*T.^2;
    k_IncoUp = func_k_Inco((dat.T_Cu2 + dat.T_Inco1)/2); % W m-2 K-1
    k_IncoDown = func_k_Inco((dat.T_Inco2 + dat.T_Cu3)/2); % W m-2 K-1
    func_k_H25 = @(T) 9.9905357 + .0205437*T -3e-6*T.^2;
    k_H25 = func_k_H25((dat.T_Inco1 + dat.T_Inco2)/2); % W m-2 K-1
    
    % delta T
    dT_IncoUp = dat.T_Cu2 - dat.T_Inco1; % K
    dT_IncoDown = dat.T_Inco2 - dat.T_Cu3; % K
    dT_H25 = dat.T_Inco1 - dat.T_Inco2; % K
    
    % heat fluxes
    t_Cu = 1.5e-3; % thickness, m
    t_Inco = 6.5e-3; % m
    h_IncoUp = ((t_Cu/k_Cu) + (t_Inco./k_IncoUp)).^-1; % htc, W m-2 K-1
    h_IncoDown = ((t_Cu/k_Cu) + (t_Inco./k_IncoDown)).^-1; % W m-2 K-1
    Q_IncoUp = h_IncoUp.*dT_IncoUp; % heat flux, W m-2
    Q_IncoDown = h_IncoDown.*dT_IncoDown; % W m-2
    Q_H25 = (Q_IncoUp + Q_IncoDown)/2; % W m-2
    
    % resistance r = 1/h
    t_Inco = 1.5e-3; % m
    t_H25 = 5.0e-3; % m
    h_cond = (t_Inco./k_IncoUp) + (t_H25./k_H25) + (t_Inco./k_IncoDown); % W m-2 K-1
    r = .5*((dT_H25./Q_H25)-h_cond); % m2 K W-1
    
    p = -dat.load*9.81/(16.25e-3^2*pi/4);
    alpha = .2;
    subplot(1, 2, 1)
    plot(p, r, '--', Color=[cmap(ff, :) alpha])
    subplot(1, 2, 2)
    plot(p, r, '--', Color=cmap(ff, :))

%     % fit data to y = ae^bx + c
%     f = @(c, x) c(1).*exp(c(2).*x) + c(3);
%     mdl = fitnlm(p, r, f, [1e-5 -1e-5 3e-4]);
%     C = mdl.Coefficients.Estimate;
%     R2 = mdl.Rsquared.Ordinary;
%     
%     p_plot = linspace(0, 2.5e5, 100);
%     plot(p_plot, f(C, p_plot), '-', Color=cmap(file, :))
%     
%     eqn = sprintf('y = %.2e*exp(%.2e x) + %.2e, R^2 = %.4f', C(1), C(2), C(3), R2);
%     txt{2*file-1} = lbl(file);
%     txt{2*file} = eqn;
end

subplot(1, 2, 1)
legend([file file], NumColumns=2)
subplot(1, 2, 2)
legend(file)
sgtitle("Solid = transient (r = 2); dashed = pseudo-steady")