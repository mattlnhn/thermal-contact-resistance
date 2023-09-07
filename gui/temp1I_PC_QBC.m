function [Ti] = temp1I_PC_QBC(qL, qR, Tavg, Ti, dt, geom, mat, A)
% Direct solver; 1 interface, perfect conduction, heat flux BCs

    % nodes
    N1 = geom{3}(1);
    N2 = geom{3}(2);
    N = N1 + N2;
    dx = geom{1};
    
    % physical properties based on average of boundary temps
    k1 = mat(1, 1, 1) + mat(1, 2, 1)*Tavg + mat(1, 3, 1)*Tavg^2;
    c_p1 = mat(2, 1, 1) + mat(2, 2, 1)*Tavg + mat(2, 3, 1)*Tavg^2;
    rho1 = mat(3, 1, 1) + mat(3, 2, 1)*Tavg + mat(3, 3, 1)*Tavg^2;
    alpha1 = k1*rho1^-1*c_p1^-1;
    tau1 = dt*alpha1*dx^-2;
    
    k2 = mat(1, 1, 2) + mat(1, 2, 2)*Tavg + mat(1, 3, 2)*Tavg^2;
    c_p2 = mat(2, 1, 2) + mat(2, 2, 2)*Tavg + mat(2, 3, 2)*Tavg^2;
    rho2 = mat(3, 1, 2) + mat(3, 2, 2)*Tavg + mat(3, 3, 2)*Tavg^2;
    alpha2 = k2*rho2^-1*c_p2^-1;
    tau2 = dt*alpha2*dx^-2;
    
    % constructing matrices
    % tridiagonal A precalculated
    % boundary condition b
    b = zeros(N, 1);
    b(1) = qL*dx/k1;
    b(end) = -qR*dx/k2;
    
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


