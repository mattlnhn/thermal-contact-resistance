function [Ti] = temp2I_CR_QBC(qL, qR, Tavg, Ti, dt, h, geom, mat, A)
% Direct solver; 2 interface, contact resistance, heat flux BC

    % rename vars
    hu = h(1);
    hd = h(2);
    
    % nodes
    N1 = geom{3}(1);
    N2 = geom{3}(2);
    N3 = geom{3}(3);
    N = N1 + N2 + N3;
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
    
    k3 = mat(1, 1, 3) + mat(1, 2, 3)*Tavg + mat(1, 3, 3)*Tavg^2;
    c_p3 = mat(2, 1, 3) + mat(2, 2, 3)*Tavg + mat(2, 3, 3)*Tavg^2;
    rho3 = mat(3, 1, 3) + mat(3, 2, 3)*Tavg + mat(3, 3, 3)*Tavg^2;
    alpha3 = k3*rho3^-1*c_p3^-1;
    tau3 = dt*alpha3*dx^-2;
    
    % build diagonal tau matrix using linear indexing
    tau = eye(N);
    tau(1:N+1:N*N1) = tau1;
    tau(N*N1+N1+1:N+1:N*(N1+N2)) = tau2;
    tau(N*(N1+N2)+N1+N2+1:N+1:N*N) = tau3;
    
    % constructing matrices
    % tridiagonal A precalculated
    % boundary condition b
    b = zeros(N, 1);
    b(1) = qL*dx/k1;
    b(end) = -qR*dx/k3;
    
    % error message if stability criterion not met
    if tau1 > .5 || tau2 > .5 || tau3 > .5
        fprintf('Error: Unstable. tau1 = %d, tau2 = %d, tau3 = %d\n', tau1, tau2, tau3)
    end
    
    % modification of matrix A on interface with contact resistance
    phi1 = 2*k1*hu^-1*dx^-1;
    phi21 = 2*k2*hu^-1*dx^-1;
    phi22 = 2*k2*hd^-1*dx^-1;
    phi3 = 2*k3*hd^-1*dx^-1;
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

