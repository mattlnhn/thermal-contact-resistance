function [Ti] = twoInterfaceCR(TL, TR, Ti, dt, H, geom, mat)
%TWOINTERFACECR Summary of this function goes here
% TL                            left BC temperature
% TR                            right BC temperature
% Ti                            temperature distribution at last time step
% dt                            time step
% H                             HTC at interface with contact resistance
% geom          geom.dx_d       distance btwn nodes in downstream section
%               geom.L1_d       length of 1st downstream section
%               geom.L2_d       length of 2nd downstream section
%               geom.L3_d       length of 3rd downstream section
% mat           mat.Fk1         anonymous function for k1(T)
%               mat.Fc_p1       anonymous function for c_p1(T)
%               mat.Frho1       anonymous function for rho1(T)
%               mat.Fk2         etc.
%               mat.Fc_p2
%               mat.Frho2
%               mat.Fk3
%               mat.Fc_p3
%               mat.Frho3

% geometry
geom.Ltotal = geom.L1_d + geom.L2_d + geom.L3_d; % total length

% nodes
N = geom.Ltotal/geom.dx_d;
N1 = geom.L1_d/geom.dx_d;
N2 = geom.L2_d/geom.dx_d;
N3 = geom.L3_d/geom.dx_d;

% constructing matrices
e = ones(N, 1);
A = spdiags([e -2*e e], [-1 0 1], N, N);
A = full(A);
A(1, 1) = -3; A(end, end) = -3;
b = zeros(N, 1);
b(1) = 2*TL;
b(end) = 2*TR;

% physical properties based on average of boundary temps
Tavg = .5*(TL + TR);

k1 = mat.Fk1(Tavg);
c_p1 = mat.Fc_p1(Tavg);
rho1 = mat.Frho1(Tavg);
alpha1 = k1*rho1^-1*c_p1^-1;
tau1 = dt*alpha1*geom.dx_d^-2;
tau1D = tau1*eye(N1); % diag matrix of tau1

k2 = mat.Fk2(Tavg);
c_p2 = mat.Fc_p2(Tavg);
rho2 = mat.Frho2(Tavg);
alpha2 = k2*rho2^-1*c_p2^-1;
tau2 = dt*alpha2*geom.dx_d^-2;
tau2D = tau2*eye(N2); % diag matrix of tau2

k3 = mat.Fk3(Tavg);
c_p3 = mat.Fc_p3(Tavg);
rho3 = mat.Frho3(Tavg);
alpha3 = k3*rho3^-1*c_p3^-1;
tau3 = dt*alpha3*geom.dx_d^-2;
tau3D = tau3*eye(N3); % diag matrix of tau3

% error message if stability criterion not met
if tau1 > .5 || tau2 > .5 || tau3 > .5
    fprintf('Error: Unstable. tau1 = %d, tau2 = %d, tau3 = %d\n', tau1, tau2, tau3)
end

% block diagonal matrix so that diffusivity corresponds to material
tau = blkdiag(tau1D, tau2D, tau3D);

% modification of matrix A on interface with contact resistance
phi1 = 2*k1*H^-1*geom.dx_d^-1;
phi2 = 2*k2*H^-1*geom.dx_d^-1;
phi3 = 2*k3*H^-1*geom.dx_d^-1;
xi1_12 = 2*phi1/(1-(1+phi1)*(1+phi2));
xi2_12 = 2*phi2/(1-(1+phi1)*(1+phi2));
xi2_23 = 2*phi2/(1-(1+phi2)*(1+phi3));
xi3_23 = 2*phi3/(1-(1+phi2)*(1+phi3));

A(N1, N1) = xi2_12-1;
A(N1, N1+1) = -xi2_12;
A(N1+1, N1) = -xi1_12;
A(N1+1, N1+1) = xi1_12-1;

A(N1+N2, N1+N2) = xi3_23-1;
A(N1+N2, N1+N2+1) = -xi3_23;
A(N1+N2+1, N1+N2) = -xi2_23;
A(N1+N2+1, N1+N2+1) = xi2_23-1;

% T at next time step
Ti = tau*(A*Ti + b) + Ti;

% error message if NaN
if any(isnan(Ti))
    fprintf('Error: NaN.\n')
end

end

