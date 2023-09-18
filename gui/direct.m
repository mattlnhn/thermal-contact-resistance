function [nextTemperature] = direct(boundaryTemperatures, ...
    averageTemperatures, currentTemperature, h, dt, geometry, materials, A)
% boundaryTemperatures  1x2 matrix, temps on L&R boundary
% averageTemperatures   1x5 matrix, average temperature in each material
% currentTemperature    Nx1 matrix, temperature distribution
% h                     heat transfer coefficient
% dt                    time step
% geometry{1}           lengths of each material
% geometry{2}           total length
% geometry{3}           dx
% geometry{4}           nodes in each section
% geometry{5}           total nodes
% geometry{6}           node numbers for T_Cu2, T_Inco1, T_Inco2, T_Cu3
% materials             3x3x5 matrix
%                       rows: k; c_p; rho
%                       cols: [a b c] where f(T) = a + bT + cT^2
%                       layers: material 1, 2, 3, 4, 5
% A                     precalculated tridiagonal A

    N = geometry{5};
    N1 = geometry{4}(1);
    N2 = geometry{4}(2);
    N3 = geometry{4}(3);
    N4 = geometry{4}(4);
    dx = geometry{3};

    % avg temps
    T1 = averageTemperatures(1);
    T2 = averageTemperatures(2);
    T3 = averageTemperatures(3);
    T4 = averageTemperatures(4);
    T5 = averageTemperatures(5);

    % physical properties based on average of boundary temps
    k1 = materials(1, 1, 1) + materials(1, 2, 1)*T1 + ...
        materials(1, 3, 1)*T1^2;
    c_p1 = materials(2, 1, 1) + materials(2, 2, 1)*T1 + ...
        materials(2, 3, 1)*T1^2;
    rho1 = materials(3, 1, 1) + materials(3, 2, 1)*T1 + ...
        materials(3, 3, 1)*T1^2;
    alpha1 = k1*rho1^-1*c_p1^-1;
    tau1 = dt*alpha1*dx^-2;
    
    k2 = materials(1, 1, 2) + materials(1, 2, 2)*T2 + ...
        materials(1, 3, 2)*T2^2;
    c_p2 = materials(2, 1, 2) + materials(2, 2, 2)*T2 + ...
        materials(2, 3, 2)*T2^2;
    rho2 = materials(3, 1, 2) + materials(3, 2, 2)*T2 + ...
        materials(3, 3, 2)*T2^2;
    alpha2 = k2*rho2^-1*c_p2^-1;
    tau2 = dt*alpha2*dx^-2;
    
    k3 = materials(1, 1, 3) + materials(1, 2, 3)*T3 + ...
        materials(1, 3, 3)*T3^2;
    c_p3 = materials(2, 1, 3) + materials(2, 2, 3)*T3 + ...
        materials(2, 3, 3)*T3^2;
    rho3 = materials(3, 1, 3) + materials(3, 2, 3)*T3 + ...
        materials(3, 3, 3)*T3^2;
    alpha3 = k3*rho3^-1*c_p3^-1;
    tau3 = dt*alpha3*dx^-2;

    k4 = materials(1, 1, 4) + materials(1, 2, 4)*T4 + ...
        materials(1, 3, 4)*T4^2;
    c_p4 = materials(2, 1, 4) + materials(2, 2, 4)*T4 + ...
        materials(2, 3, 4)*T4^2;
    rho4 = materials(3, 1, 4) + materials(3, 2, 4)*T4 + ...
        materials(3, 3, 4)*T4^2;
    alpha4 = k4*rho4^-1*c_p4^-1;
    tau4 = dt*alpha4*dx^-2;

    k5 = materials(1, 1, 5) + materials(1, 2, 5)*T5 + ...
        materials(1, 3, 5)*T5^2;
    c_p5 = materials(2, 1, 5) + materials(2, 2, 5)*T5 + ...
        materials(2, 3, 5)*T5^2;
    rho5 = materials(3, 1, 5) + materials(3, 2, 5)*T5 + ...
        materials(3, 3, 5)*T5^2;
    alpha5 = k5*rho5^-1*c_p5^-1;
    tau5 = dt*alpha5*dx^-2;

    % build diagonal tau
    tau = eye(N);
    tau(1:N+1:N*N1) = tau1;
    tau(N*N1+N1+1:N+1:N*(N1+N2)) = tau2;
    tau(N*(N1+N2)+(N1+N2)+1:N+1:N*N) = tau3;
    tau(N*(N1+N2+N3)+(N1+N2+N3)+1:N+1:N*N) = tau4;
    tau(N*(N1+N2+N3+N4)+(N1+N2+N3+N4)+1:N+1:N*N) = tau5;

    %fprintf("tau1 %d tau2 %d tau3 %d tau4 %d tau5 %d\n", tau1, tau2, tau3, tau4, tau5)

    % constructing matrices
    % tridiagonal A precalculated
    % boundary condition b
    b = zeros(N, 1);
    %b(1) = qL*dx/k1;
    %b(end) = -qR*dx/k3;
    b(1) = 2*boundaryTemperatures(1);
    b(end) = 2*boundaryTemperatures(end);

    % boundary 1->2, assumed perfect conduction
    kappa12 = k1/k2;
    gamma12 = (kappa12-1)/(kappa12+1);
    A(N1, N1) = -(2-gamma12);
    A(N1, N1+1) = 2/(kappa12+1);
    A(N1+1, N1) = 2*kappa12/(kappa12+1);
    A(N1+1, N1+1) = -(2+gamma12);

    % boundary 2->3, contact resistance
    phi2_23 = 2*k2*h^-1*dx^-1;
    phi3_23 = 2*k3*h^-1*dx^-1;
    xi2_23 = 2*phi2_23/(1-(1+phi2_23)*(1+phi3_23));
    xi3_23 = 2*phi3_23/(1-(1+phi2_23)*(1+phi3_23));
    A(N1+N2, N1+N2) = xi3_23-1;
    A(N1+N2, N1+N2+1) = -xi3_23;
    A(N1+N2+1, N1+N2) = -xi2_23;
    A(N1+N2+1, N1+N2+1) = xi2_23-1;

    % boundary 3->4, contact resistance
    phi3_34 = 2*k3*h^-1*dx^-1;
    phi4_34 = 2*k4*h^-1*dx^-1;
    xi3_34 = 2*phi3_34/(1-(1+phi3_34)*(1+phi4_34));
    xi4_34 = 2*phi4_34/(1-(1+phi3_34)*(1+phi4_34));
    A(N1+N2+N3, N1+N2+N3) = xi4_34-1;
    A(N1+N2+N3, N1+N2+N3+1) = -xi4_34;
    A(N1+N2+N3+1, N1+N2+N3) = -xi3_34;
    A(N1+N2+N3+1, N1+N2+N3+1) = xi3_34-1;

    % boundary 4->5, assumed perfect conduction
    kappa45 = k4/k5;
    gamma45 = (kappa45-1)/(kappa45+1);
    A(N1+N2+N3+N4, N1+N2+N3+N4) = -(2-gamma45);
    A(N1+N2+N3+N4, N1+N2+N3+N4+1) = 2/(kappa45+1);
    A(N1+N2+N3+N4+1, N1+N2+N3+N4) = 2*kappa45/(kappa45+1);
    A(N1+N2+N3+N4+1, N1+N2+N3+N4+1) = -(2+gamma45);

    % update
    nextTemperature = tau*(A*currentTemperature + b) + currentTemperature;

    % error message if NaN
    if any(isnan(nextTemperature))
        fprintf('Error: NaN.\n')
    end

end

