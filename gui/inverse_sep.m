function [h] = inverse_sep(app, filename, geometry, materials, parameters)
%INVERSE Summary of this function goes here
%   Detailed explanation goes here
    %% progress bar
    fig = app.ihcpUIFigure;
    d = uiprogressdlg(fig, 'Title', 'In progress', 'Message', ...
        'Computing heat transfer coefficients...');
    
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

    initialT = interp1([1; geometry{6}'; N], ...
        [dat.T_Cu1(1); dat.T_Cu2(1); dat.T_Inco1(1); ...
        dat.T_Inco2(1); dat.T_Cu3(1); dat.T_Cu4(1)], ...
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
        Y = [dat.T_Cu2(m+1:m+r)';
            dat.T_Inco1(m+1:m+r)';
            dat.T_Inco2(m+1:m+r)';
            dat.T_Cu3(m+1:m+r)'];
        % for average temp
        y = [dat.T_Cu1(m+1:m+r)'; Y; dat.T_Cu4(m+1:m+r)'];
        Tavg = .5*(y(2:end, :)-y(1:end-1, :));
    
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
            dhu = sum((Y-Ts).*Xu, 'all')/sum(Xu.^2, 'all');
            dhd = sum((Y-Ts).*Xd, 'all')/sum(Xd.^2, 'all');
    
            h = h + [dhu dhd];
    
            erroru = sum(.5*(Y-Tdhus).^2, 'all');
            errord = sum(.5*(Y-Tdhds).^2, 'all');
    
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
        for n = 1:steps(2)+1
            initialT = direct_sep([dat.T_Cu1(m+1) dat.T_Cu4(m+1)], ...
                    Tavg(:, 1), initialT, h, dt, geometry, materials, A);
        end
    
        d.Value = m/(totalSteps-r);
    end

    %% write file
    [path, file, ext] = fileparts(filename);
    writematrix(hStore, path+"\h_"+file+ext, "WriteMode", "overwrite");

    %% tidy up
    close(d) % kill progress bar

    if app.plotCheck.Value % basic plot if asked
        figure()
        hold on
        plot(-dat.load(1:end-r)/(.25*pi*16.25e-3^2), hStore(2:end, 1).^-1, 'r')
        xlabel("Pressure [Pa]")
        ylabel("Contact resistance h^{-1} [m^2 K W^{-1}]")
        ylim([0 1.5*max(hStore.^-1, [], "all")])
        grid on
        grid minor
    end
end

