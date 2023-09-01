function [Ti] = temp2I_CR_QBC(qL, qR, Tavg, Ti, dt, h, geom, mat)
% temp distribution for 2 interface, contact resistance, heat flux BC
% QL                            left BC heat flux
% QR                            right BC heat flux
% Ti                            temperature distribution at last time step
% dt                            time step
% h                             [hu hd] htcs at upstream/downstream
% geom          geom.dx         distance btwn nodes
%               geom.L1         length of 1st section
%               geom.L2         length of 2nd section
%               geom.L3         length of 3rd section
% mat           mat.Fk1         anonymous function for k1(T)
%               mat.Fc_p1       anonymous function for c_p1(T)
%               mat.Frho1       anonymous function for rho1(T)
%               mat.Fk2         etc.
%               mat.Fc_p2
%               mat.Frho2
%               mat.Fk3
%               mat.Fc_p3
%               mat.Frho3

hu = h(1);
hd = h(2);

% geometry
geom.Ltotal = geom.L1 + geom.L2 + geom.L3; % total length

% nodes
N = geom.Ltotal/geom.dx;
N1 = geom.L1/geom.dx;
N2 = geom.L2/geom.dx;
N3 = geom.L3/geom.dx;

% physical properties based on average of boundary temps
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

k3 = mat.Fk3(Tavg);
c_p3 = mat.Fc_p3(Tavg);
rho3 = mat.Frho3(Tavg);
alpha3 = k3*rho3^-1*c_p3^-1;
tau3 = dt*alpha3*geom.dx^-2;

% build diagonal tau matrix using linear indexing
tau = eye(N);
tau(1:N+1:N*N1) = tau1;
tau(N*N1+N1+1:N+1:N*(N1+N2)) = tau2;
tau(N*(N1+N2)+N1+N2+1:N+1:N*N) = tau3;

% constructing matrices
% tridiagonal A
A = diag(ones(N-1, 1), -1) + diag(-2*ones(N, 1), 0) + diag(ones(N-1, 1), 1);
A(1, 1) = -1; A(end, end) = -1;
% boundary condition b
b = zeros(N, 1);
b(1) = qL*geom.dx/k1;
b(end) = -qR*geom.dx/k3;

% error message if stability criterion not met
if tau1 > .5 || tau2 > .5 || tau3 > .5
    fprintf('Error: Unstable. tau1 = %d, tau2 = %d, tau3 = %d\n', tau1, tau2, tau3)
end

% modification of matrix A on interface with contact resistance
phi1 = 2*k1*hu^-1*geom.dx^-1;
phi21 = 2*k2*hu^-1*geom.dx^-1;
phi22 = 2*k2*hd^-1*geom.dx^-1;
phi3 = 2*k3*hd^-1*geom.dx^-1;
xi1_12 = 2*phi1/(1-(1+phi1)*(1+phi21));
xi2_12 = 2*phi21/(1-(1+phi1)*(1+phi21));
xi2_23 = 2*phi22/(1-(1+phi22)*(1+phi3));
xi3_23 = 2*phi3/(1-(1+phi22)*(1+phi3));

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

