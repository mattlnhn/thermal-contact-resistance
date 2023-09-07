function upstreamdownstream(app, filename, geom, mat, param, dt)
% Inverse solver for heat flux in up/downstream sections
    %% progress bar
    fig = app.ihcpUIFigure;
    d = uiprogressdlg(fig, 'Title', 'In progress (1 of 3)', 'Message', ...
        'Computing heat flux in upstream section...');

    %% section geometry
    gu = geom.up;
    Nu = round(sum(gu{3}, "all"));
    gd = geom.down;
    Nd = round(sum(gd{3}, "all"));
    
    %% precalculate A matrices
    Au = diag(ones(Nu-1, 1), -1) + diag(-2*ones(Nu, 1), 0) + ...
        diag(ones(Nu-1, 1), 1);
    Au(1, 1) = -1; Au(end, end) = -1;
    
    Ad = diag(ones(Nd-1, 1), -1) + diag(-2*ones(Nd, 1), 0) + ...
        diag(ones(Nd-1, 1), 1);
    Ad(1, 1) = -1; Ad(end, end) = -1;
    
    %% material properties
    mu = mat.up;
    md = mat.down;
    
    %% parameters
    r = param.r;
    epsilon = 1e-2;
    RTOLq = param.RTOLqh;
    RTOLerror = param.RTOLerror;
    maxIter = param.maxiter;
    qInitial = param.qInitial;
    
    %% import
    dat = readtable(filename);
    totalSteps = length(dat.time);
    qu = zeros(totalSteps-r-1, 2);
    qd = qu;
    qu(1, :) = qInitial;
    qd(1, :) = qInitial;
    
    %% time iteration upstream
    % initial temp distribution
    initialT = interp1([1; Nu], [dat.T_Cu2(1); dat.T_Inco1(1)], (1:Nu)');
    
    for m = 1:totalSteps-r-1
        % reset counters
        converged = 0;
        iter = 1;
        errorL = 1e20;
        errorR = errorL;
    
        % time data
        time = dat.time(m:m+r);
        % no. of steps of dt from initial time
        steps = round((time - time(1))/dt);
    
        % temp matrices incl. initial time
        T = zeros(Nu, steps(end)+1);
        T(:, 1) = initialT;
        TdqL = T;
        TdqR = T;
    
        % start with previous values of q
        qL = qu(m, 1);
        qR = qu(m, 2);
    
        % measured data r steps into future
        Y = [dat.T_Cu2(m:m+r)';
            dat.T_Inco1(m:m+r)'];
        Tavg = mean(Y, 1);
    
        while converged == 0
            preverrorL = errorL;
            preverrorR = errorR;
    
            for n = 1:steps(end)
                % index for avg temp
                ai = sum(steps<n);
                % temp updates
                T(:, n+1) = temp1I_PC_QBC(qL, qR, Tavg(ai), T(:, n), dt, ...
                    gu, mu, Au);
                TdqL(:, n+1) = temp1I_PC_QBC((1+epsilon)*qL, qR, Tavg(ai), ...
                    TdqL(:, n), dt, gu, mu, Au);
                TdqR(:, n+1) = temp1I_PC_QBC(qL, (1+epsilon)*qR, Tavg(ai), ...
                    TdqR(:, n), dt, gu, mu, Au);
            end
    
            % T at sensor locations & measurement times
            Ts = [T(1, steps+1); T(end, steps+1)];
            TdqLs = [TdqL(1, steps+1); TdqL(end, steps+1)];
            TdqRs = [TdqR(1, steps+1); TdqR(end, steps+1)];
    
            % sensitivity coefficients
            XL = (TdqLs - Ts)./(epsilon*qL);
            XR = (TdqRs - Ts)./(epsilon*qR);
    
            % newton method increment
            dqL = sum((Y-Ts).*XL, 'all')/sum(XL.^2, 'all');
            dqR = sum((Y-Ts).*XR, 'all')/sum(XR.^2, 'all');
    
            qL = qL + dqL;
            qR = qR + dqR;
    
            errorL = sum(.5*(Y-TdqLs).^2, 'all');
            errorR = sum(.5*(Y-TdqRs).^2, 'all');
    
            % convergence criteria
            % relative change in q
            if abs(dqL/qL) < RTOLq && abs(dqR/qR) < RTOLq
                converged = 1;
            end
            % relative change in error
            if abs((errorL-preverrorL)/preverrorL) < RTOLerror && ...
                    abs((errorR-preverrorR)/preverrorR) < RTOLerror
                converged = 1;
            end
            % timeout
            if iter > maxIter
                converged = 1;
            end
    
            iter = iter + 1;
    
        end
    
        % save q
        qu(m+1, :) = [qL qR];
        
        % next step initial temp
        for n = 1:steps(2)
            initialT = temp1I_PC_QBC(qL, qR, Tavg(1), initialT, dt, gu, mu, Au);
        end
        
        d.Value = m/(totalSteps-r-1);
    end
    
    %% time iteration downstream
    % initial temp distribution
    initialT = interp1([1; Nd], [dat.T_Inco2(1); dat.T_Cu3(1)], (1:Nd)');
    
    % rename progress bar
    d = uiprogressdlg(fig, 'Title', 'In progress (2 of 3)', 'Message', ...
        'Computing heat flux in downstream section...');

    for m = 1:totalSteps-r-1
        % reset counters
        converged = 0;
        iter = 1;
        errorL = 1e20;
        errorR = errorL;
    
        % time data
        time = dat.time(m:m+r);
        steps = round((time - time(1))/dt);
    
        % temp matrices incl. initial time
        T = zeros(Nd, steps(end)+1);
        T(:, 1) = initialT;
        TdqL = T;
        TdqR = T;
    
        % start with previous values of q
        qL = qd(m, 1);
        qR = qd(m, 2);
    
        % measured data r steps into future
        Y = [dat.T_Inco2(m:m+r)';
            dat.T_Cu3(m:m+r)'];
        Tavg = mean(Y, 1);
    
        while converged == 0
            preverrorL = errorL;
            preverrorR = errorR;
    
            for n = 1:steps(end)
                % index for avg temp
                ai = sum(steps<n);
                % temp updates
                T(:, n+1) = temp1I_PC_QBC(qL, qR, Tavg(ai), T(:, n), dt, ...
                    gd, md, Ad);
                TdqL(:, n+1) = temp1I_PC_QBC((1+epsilon)*qL, qR, Tavg(ai), ...
                    TdqL(:, n), dt, gd, md, Ad);
                TdqR(:, n+1) = temp1I_PC_QBC(qL, (1+epsilon)*qR, Tavg(ai), ...
                    TdqR(:, n), dt, gd, md, Ad);
            end
    
            % T at sensor locations & measurement times
            Ts = [T(1, steps+1); T(end, steps+1)];
            TdqLs = [TdqL(1, steps+1); TdqL(end, steps+1)];
            TdqRs = [TdqR(1, steps+1); TdqR(end, steps+1)];
    
            % sensitivity coefficients
            XL = (TdqLs - Ts)./(epsilon*qL);
            XR = (TdqRs - Ts)./(epsilon*qR);
    
            % newton method increment
            dqL = sum((Y-Ts).*XL, 'all')/sum(XL.^2, 'all');
            dqR = sum((Y-Ts).*XR, 'all')/sum(XR.^2, 'all');
    
            qL = qL + dqL;
            qR = qR + dqR;
    
            errorL = sum(.5*(Y-TdqLs).^2, 'all');
            errorR = sum(.5*(Y-TdqRs).^2, 'all');
    
            % convergence criteria
            % relative change in q
            if abs(dqL/qL) < RTOLq && abs(dqR/qR) < RTOLq
                converged = 1;
            end
            % relative change in error
            if abs((errorL-preverrorL)/preverrorL) < RTOLerror && ...
                    abs((errorR-preverrorR)/preverrorR) < RTOLerror
                converged = 1;
            end
            % timeout
            if iter > maxIter
                converged = 1;
            end
    
            iter = iter + 1;
    
        end
    
        % save q
        qd(m+1, :) = [qL qR];
        
        % next step initial temp
        for n = 1:steps(2)
            initialT = temp1I_PC_QBC(qL, qR, Tavg(1), initialT, dt, gd, md, Ad);
        end
    
        d.Value = m/(totalSteps-r-1);
    end
    
    % save relevant boundaries
    qu = qu(:, 2);
    qd = qd(:, 1);
    newfilename = filename+".mat";
    save(newfilename, "qu", "qd")

end

