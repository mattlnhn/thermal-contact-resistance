function [Ti] = oneInterfaceNoCR(TL, TR, Ti, dt, geom, mat)
%ONEINTERFACE Summary of this function goes here
% TL                            left BC temperature
% TR                            right BC temperature
% Ti                            temperature distribution at last time step
% dt                            time step
% geom          geom.dx_u       distance btwn nodes in upstream section
%               geom.L1_u       length of 1st upstream section
%               geom.L2_u       length of 2nd upstream section
% mat           mat.Fk0         anonymous function for k0(T)
%               mat.Fc_p0       anonymous function for c_p0(T)
%               mat.Frho0       anonymous function for rho0(T)
%               mat.Fk1         etc.
%               mat.Fc_p1
%               mat.Frho1

% geometry
geom.Ltotal = geom.L1_u + geom.L2_u; % total section length

% nodes
N = geom.Ltotal/geom.dx_u;
N0 = geom.L1_u/geom.dx_u;
N1 = geom.L2_u/geom.dx_u;

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

k0 = mat.Fk0(Tavg);
c_p0 = mat.Fc_p0(Tavg);
rho0 = mat.Frho0(Tavg);
alpha0 = k0*rho0^-1*c_p0^-1;
tau0 = dt*alpha0*geom.dx_u^-2;
tau0D = tau0*eye(N0); % diag matrix of tau1

k1 = mat.Fk1(Tavg);
c_p1 = mat.Fc_p1(Tavg);
rho1 = mat.Frho1(Tavg);
alpha1 = k1*rho1^-1*c_p1^-1;
tau1 = dt*alpha1*geom.dx_u^-2;
tau1D = tau1*eye(N1); % diag matrix of tau1

% error message if stability criterion not met
if tau0 > .5 || tau1 > .5
    fprintf('Error: Unstable. tau0 = %d, tau1 = %d\n', tau0, tau1)
end

% block diagonal matrix so that diffusivity corresponds to material
tau = blkdiag(tau0D, tau1D);

% modification of matrix A on interface with no contact resistance
kappa = k0/k1;
gamma = (kappa-1)/(kappa+1);
A(N0, N0) = -(2-gamma); A(N0+1, N0+1) = -(2+gamma);
A(N0, N0+1) = 2/(kappa+1); A(N0+1, N0) = 2*kappa/(kappa+1);

% T at next time step
Ti = tau*(A*Ti + b) + Ti;

% error message if NaN
if any(isnan(Ti))
    fprintf('Error: NaN.\n')
end

end

