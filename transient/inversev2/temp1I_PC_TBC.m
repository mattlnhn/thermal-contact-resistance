function [Ti] = temp1I_PC_TBC(TL, TR, Ti, dt, geom, mat)
% temp distribution for 1 interface, perfect conduction, temperature BCs
% TL                            left BC temperature
% TR                            right BC temperature
% Ti                            temperature distribution at last time step
% dt                            time step
% geom          geom.dx         distance btwn nodes in upstream section
%               geom.L1         length of 1st section
%               geom.L2         length of 2nd section
% mat           mat.Fk1         anonymous function for k1(T)
%               mat.Fc_p1       anonymous function for c_p1(T)
%               mat.Frho1       anonymous function for rho1(T)
%               mat.Fk1         etc.
%               mat.Fc_p1
%               mat.Frho1

% geometry
geom.Ltotal = geom.L1 + geom.L2; % total section length

% nodes
N = geom.Ltotal/geom.dx;
N1 = geom.L1/geom.dx;
N2 = geom.L2/geom.dx;

% constructing matrices
% tridiagonal A
A = diag(ones(N-1, 1), -1) + diag(-2*ones(N, 1), 0) + diag(ones(N-1, 1), 1);
A(1, 1) = -3; A(end, end) = -3;
% boundary condition b
b = zeros(N, 1);
b(1) = 2*TL;
b(end) = 2*TR;

% physical properties based on average of boundary temps
Tavg = .5*(TL + TR);

k1 = mat.Fk1(Tavg);
c_p1 = mat.Fc_p1(Tavg);
rho1 = mat.Frho1(Tavg);
alpha1 = k1*rho1^-1*c_p1^-1;
tau1 = dt*alpha1*geom.dx^-2;

k2 = mat.Fk2(Tavg);
c_p2 = mat.Fc_p2(Tavg);
rho2 = mat.Frho2(Tavg);
alpha2 = k2*rho2^-1*c_p2^-1;
tau2 = dt*alpha2*geom.dx^-2;

% build diagonal tau matrix using linear indexing
tau = eye(N);
tau(1:N+1:N*N1) = tau1;
tau(N*N1+N1+1:N+1:N*N) = tau2;

% error message if stability criterion not met
if tau1 > .5 || tau2 > .5
    fprintf('Error: Unstable. tau1 = %d, tau2 = %d\n', tau1, tau2)
end

% modification of matrix A on interface with no contact resistance
kappa = k1/k2;
gamma = (kappa-1)/(kappa+1);
A(N1, N1) = -(2-gamma); A(N1+1, N1+1) = -(2+gamma);
A(N1, N1+1) = 2/(kappa+1); A(N1+1, N1) = 2*kappa/(kappa+1);

% T at next time step
Ti = tau*(A*Ti + b) + Ti;

% error message if NaN
if any(isnan(Ti))
    fprintf('Error: NaN.\n')
end

end

