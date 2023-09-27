function [h] = inverse_sepextra(filename, geometry, materials, parameters)
%INVERSE Summary of this function goes here
%   Detailed explanation goes here
    
    %% precalculate A matrix
    N = geometry{5};
    A = diag(ones(N-1, 1), -1) + diag(-2*ones(N, 1), 0) + ...
        diag(ones(N-1, 1), 1);
    A(1, 1) = -3; A(end, end) = -3;

    %% parameters
    r = parameters.r;
    dt = parameters.dt;
    epsilon = parameters.epsilon;
    RTOLh = parameters.RTOLh;
    RTOLerror = parameters.RTOLerror;
    maxIter = parameters.maxiter;
    hInitial = parameters.hInitial;

    %% import
    dat = readtable(filename);
    totalSteps = length(dat.time);
    hStore = zeros(totalSteps-r, 2);
    hStore(1, :) = hInitial;

%     initialT = interp1([1; geometry{6}'; N], ...
%         [dat.T_Cu1(1); dat.T_Cu2(1); dat.T_Inco1a(1); dat.T_Inco1(1); ...
%         dat.T_H25a(1); dat.T_H25b(1); dat.T_Inco2(1); dat.T_Inco2a(1); ...
%         dat.T_Cu3(1); dat.T_Cu4(1)], ...
%         (1:N)', "linear", "extrap");
    initialT = interp1([1; geometry{6}'; N], ...
        [dat.T_Cu1(1); dat.T_Cu2(1); dat.T_Inco1(1); ...
        dat.T_H25(1); dat.T_Inco2(1); ...
        dat.T_Cu3(1); dat.T_Cu4(1)], ...
        (1:N)', "linear", "extrap");

    for m = 1:totalSteps-r
        % reset counters
        converged = 0;
        iter = 1;
        erroru = 1e20;
        errord = erroru;
    
        % time data
        time = dat.time(m:m+r);
        % no. of steps of dt from initial time
        steps = round((time - time(1))/dt);
    
        % temp matrices incl. initial time
        T = zeros(N, steps(end)+1);
        T(:, 1) = initialT;
        Tdhu = T;
        Tdhd = T;
    
        % start with previous values of h
        h = hStore(m, :);
    
        % measured data r steps into future
%         Y = [dat.T_Cu2(m+1:m+r)';
%             dat.T_Inco1a(m+1:m+r)';
%             dat.T_Inco1(m+1:m+r)';
%             dat.T_H25a(m+1:m+r)';
%             dat.T_H25b(m+1:m+r)';
%             dat.T_Inco2(m+1:m+r)';
%             dat.T_Inco2a(m+1:m+r)';
%             dat.T_Cu3(m+1:m+r)'];
        Y = [dat.T_Cu2(m+1:m+r)';
            dat.T_Inco1(m+1:m+r)';
            dat.T_H25(m+1:m+r)';
            dat.T_Inco2(m+1:m+r)';
            dat.T_Cu3(m+1:m+r)'];
        % for average temp
        y = [dat.T_Cu1(m+1:m+r)'; Y; dat.T_Cu4(m+1:m+r)'];
        Tavg = .5*(y(2:end, :)+y(1:end-1, :));
    
        while converged == 0
            preverroru = erroru;
            preverrord = errord;
    
            for n = 1:steps(end)
                % index for avg temp
                ai = sum(steps<n);
                % temp updates
                T(:, n+1) = direct_sep([dat.T_Cu1(m+ai) dat.T_Cu4(m+ai)], ...
                    Tavg(:, ai), T(:, n), h, dt, geometry, materials, A);
                Tdhu(:, n+1) = direct_sep([dat.T_Cu1(m+ai) dat.T_Cu4(m+ai)], ...
                    Tavg(:, ai), Tdhu(:, n), [(1+epsilon)*h(1) h(2)], dt, geometry, ...
                    materials, A);
                Tdhd(:, n+1) = direct_sep([dat.T_Cu1(m+ai) dat.T_Cu4(m+ai)], ...
                    Tavg(:, ai), Tdhd(:, n), [h(1) (1+epsilon)*h(2)], dt, geometry, ...
                    materials, A);
            end
    
            % T at sensor locations & measurement times, excl. initial
            Ts = [T(geometry{6}, steps(2:end)+1)];
            Tdhus = [Tdhu(geometry{6}, steps(2:end)+1)];
            Tdhds = [Tdhd(geometry{6}, steps(2:end)+1)];

            % sensitivity coefficients
            Xu = (Tdhus - Ts)./(epsilon*h(1));
            Xd = (Tdhds - Ts)./(epsilon*h(2));
    
            % newton method increment
            dhu = sum((Y-Tdhus).*Xu, 'all')/sum(Xu.^2, 'all');
            dhd = sum((Y-Tdhds).*Xd, 'all')/sum(Xd.^2, 'all');
    
            h = h + [dhu dhd];
    
            erroru = sum(.5*(Y-Ts).^2, 'all');
            errord = sum(.5*(Y-Ts).^2, 'all');
    
            % convergence criteria
            % relative change in h
            if abs(dhu/h(1)) < RTOLh && abs(dhd/h(2)) < RTOLh
                converged = 1;
            end
            % relative change in error
            if abs((erroru-preverroru)/preverroru) < RTOLerror && ...
                    abs((errord-preverrord)/preverrord) < RTOLerror
                converged = 1;
            end
            % timeout
            if iter > maxIter
                converged = 1;
            end
    
            iter = iter + 1;
    
        end
    
        % save q
        hStore(m+1, :) = h;
        
        % next step initial temp
        for n = 1:steps(2)
            initialT = direct_sep([dat.T_Cu1(m+1) dat.T_Cu4(m+1)], ...
                    Tavg(:, 1), initialT, h, dt, geometry, materials, A);
        end
    
        fprintf("%.2f complete\n", 100*m/(totalSteps-r));
    end

    %% write file
    [path, file, ext] = fileparts(filename);
    writematrix(hStore, path+"\h_"+file+ext, "WriteMode", "overwrite");

end

