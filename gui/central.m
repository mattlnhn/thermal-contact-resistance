function central(app, filename, geom, mat, param, dt)
% Inverse solver for unknown heat transfer coefficients at boundaries
    %% progress bar
    fig = app.ihcpUIFigure;
    d = uiprogressdlg(fig, 'Title', 'In progress (3 of 3)', 'Message', ...
        'Computing heat transfer coefficients in central section...');

    %% section geometry
    gc = geom.cent;
    Nc = round(sum(gc{3}, "all"));
    TCH25 = floor(Nc*(gc{2}(1)+gc{4})/sum(gc{2}, "all"));
    
    %% precalculate A matrices
    Ac = diag(ones(Nc-1, 1), -1) + diag(-2*ones(Nc, 1), 0) + ...
        diag(ones(Nc-1, 1), 1);
    Ac(1, 1) = -1; Ac(end, end) = -1;
    
    %% material properties
    mc = mat.cent;
    
    %% precaulculated upstream/downstream data
    load(filename+".mat", "qu", "qd")
  
    %% parameters
    r = param.r;
    epsilon = 1e-6;
    RTOLh = param.RTOLqh;
    RTOLerror = param.RTOLerror;
    maxIter = param.maxiter;
    hInitial = param.hInitial;

    %% import
    dat = readtable(filename);
    totalSteps = length(dat.time);
    hStore = zeros(totalSteps-r-1, 2);
    hStore(1, :) = hInitial;
    
    %% time iteration
    % initial temp distribution
    initialT = interp1([1; TCH25; Nc], [dat.T_Inco1(1); dat.T_H25(1); dat.T_Inco2(1)], (1:Nc)');
    
    if app.IndividualButton.Value == 1 % individual htcs

        for m = 1:totalSteps-r-1
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
            T = zeros(Nc, steps(end)+1);
            T(:, 1) = initialT;
            Tdhu = T;
            Tdhd = T;
        
            % start with previous values of h
            hu = hStore(m, 1);
            hd = hStore(m, 2);
        
            % measured data r steps into future
            Y = [dat.T_Inco1(m:m+r)';
                dat.T_H25(m:m+r)';
                dat.T_Inco2(m:m+r)'];
            Tavg = mean(Y, 1);
        
            while converged == 0
                preverroru = erroru;
                preverrord = errord;
        
                for n = 1:steps(end)
                    % index for avg temp
                    ai = sum(steps<n);
                    % temp updates
                    T(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                        T(:, n), dt, [hu hd], gc, mc, Ac);
                    Tdhu(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                        Tdhu(:, n), dt, [(1+epsilon)*hu hd], gc, mc, Ac);
                    Tdhd(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                        Tdhd(:, n), dt, [hu (1+epsilon)*hd], gc, mc, Ac);
                end
        
                % T at sensor locations & measurement times, excl. initial
                Ts = [T(1, steps+1); T(TCH25, steps+1); T(end, steps+1)];
                Tdhus = [Tdhu(1, steps+1); Tdhu(TCH25, steps+1); Tdhu(end, steps+1)];
                Tdhds = [Tdhd(1, steps+1); Tdhd(TCH25, steps+1); Tdhd(end, steps+1)];
        
                % sensitivity coefficients
                Xu = (Tdhus - Ts)./(epsilon*hu);
                Xd = (Tdhds - Ts)./(epsilon*hd);
        
                % newton method increment
                dhu = sum((Y-Ts).*Xu, 'all')/sum(Xu.^2, 'all');
                dhd = sum((Y-Ts).*Xd, 'all')/sum(Xd.^2, 'all');
        
                hu = hu + dhu;
                hd = hd + dhd;
        
                erroru = sum(.5*(Y-Tdhus).^2, 'all');
                errord = sum(.5*(Y-Tdhds).^2, 'all');
        
                % convergence criteria
                % relative change in h
                if abs(dhu/hu) < RTOLh && abs(dhd/hd) < RTOLh
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
            hStore(m+1, 1) = hu;
            hStore(m+1, 2) = hd;
            
            % next step initial temp
            for n = 1:steps(2)
                initialT = temp2I_CR_QBC(qu(m), qd(m), Tavg(1), ...
                        initialT, dt, [hu hd], gc, mc, Ac);
            end
        
            d.Value = m/(totalSteps-r-1);
        end

    elseif app.centreTC == 1 % avg htc, including centre tc

        for m = 1:totalSteps-r-1
            % reset counters
            converged = 0;
            iter = 1;
            error = 1e20;
        
            % time data
            time = dat.time(m:m+r);
            % no. of steps of dt from initial time
            steps = round((time - time(1))/dt);
        
            % temp matrices incl. initial time
            T = zeros(Nc, steps(end)+1);
            T(:, 1) = initialT;
            Tdh = T;
        
            % start with previous values of h
            h = hStore(m, 1);
        
            % measured data r steps into future
            Y = [dat.T_Inco1(m:m+r)';
                dat.T_H25(m:m+r)';
                dat.T_Inco2(m:m+r)'];
            Tavg = mean(Y, 1);
        
            while converged == 0
                preverror = error;
        
                for n = 1:steps(end)
                    % index for avg temp
                    ai = sum(steps<n);
                    % temp updates
                    T(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                        T(:, n), dt, [h h], gc, mc, Ac);
                    Tdh(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                        Tdh(:, n), dt, (1+epsilon)*[h h], gc, mc, Ac);
                end
        
                % T at sensor locations & measurement times, excl. initial
                Ts = [T(1, steps+1); T(TCH25, steps+1); T(end, steps+1)];
                Tdhs = [Tdh(1, steps+1); Tdh(TCH25, steps+1); Tdh(end, steps+1)];

                % sensitivity coefficients
                X = (Tdhs - Ts)./(epsilon*h);
        
                % newton method increment
                dh = sum((Y-Ts).*X, 'all')/sum(X.^2, 'all');
        
                h = h + dh;
        
                error = sum(.5*(Y-Tdhs).^2, 'all');
        
                % convergence criteria
                % relative change in h
                if abs(dh/h) < RTOLh
                    converged = 1;
                end
                % relative change in error
                if abs((error-preverror)/preverror) < RTOLerror
                    converged = 1;
                end
                % timeout
                if iter > maxIter
                    converged = 1;
                end
        
                iter = iter + 1;
        
            end
        
            % save q
            hStore(m+1, 1) = h;
            hStore(m+1, 2) = h;
            
            % next step initial temp
            for n = 1:steps(2)
                initialT = temp2I_CR_QBC(qu(m), qd(m), Tavg(1), ...
                        initialT, dt, [h h], gc, mc, Ac);
            end
        
            d.Value = m/(totalSteps-r-1);
        end

    else % avg htc, excluding centre tc

        for m = 1:totalSteps-r-1
            % reset counters
            converged = 0;
            iter = 1;
            error = 1e20;
        
            % time data
            time = dat.time(m:m+r);
            % no. of steps of dt from initial time
            steps = round((time - time(1))/dt);
        
            % temp matrices incl. initial time
            T = zeros(Nc, steps(end)+1);
            T(:, 1) = initialT;
            Tdh = T;
        
            % start with previous values of h
            h = hStore(m, 1);
        
            % measured data r steps into future
            Y = [dat.T_Inco1(m:m+r)';
                dat.T_Inco2(m:m+r)'];
            Tavg = mean(Y, 1);
        
            while converged == 0
                preverror = error;
        
                for n = 1:steps(end)
                    % index for avg temp
                    ai = sum(steps<n);
                    % temp updates
                    T(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                        T(:, n), dt, [h h], gc, mc, Ac);
                    Tdh(:, n+1) = temp2I_CR_QBC(qu(m), qd(m), Tavg(ai), ...
                        Tdh(:, n), dt, (1+epsilon)*[h h], gc, mc, Ac);
                end
        
                % T at sensor locations & measurement times, excl. initial
                Ts = [T(1, steps+1); T(end, steps+1)];
                Tdhs = [Tdh(1, steps+1); Tdh(end, steps+1)];

                % sensitivity coefficients
                X = (Tdhs - Ts)./(epsilon*h);
        
                % newton method increment
                dh = sum((Y-Ts).*X, 'all')/sum(X.^2, 'all');
        
                h = h + dh;
        
                error = sum(.5*(Y-Tdhs).^2, 'all');
        
                % convergence criteria
                % relative change in h
                if abs(dh/h) < RTOLh
                    converged = 1;
                end
                % relative change in error
                if abs((error-preverror)/preverror) < RTOLerror
                    converged = 1;
                end
                % timeout
                if iter > maxIter
                    converged = 1;
                end
        
                iter = iter + 1;
        
            end
        
            % save q
            hStore(m+1, 1) = h;
            hStore(m+1, 2) = h;
            
            % next step initial temp
            for n = 1:steps(2)
                initialT = temp2I_CR_QBC(qu(m), qd(m), Tavg(1), ...
                        initialT, dt, [h h], gc, mc, Ac);
            end
        
            d.Value = m/(totalSteps-r-1);
        end

    end
    
    %% write file
    [path, file, ext] = fileparts(filename);
    writematrix(hStore, path+"\h_"+file+ext, "WriteMode", "overwrite");

    %% tidy up
    close(d) % kill progress bar

    if app.plotCheck.Value % basic plot if asked
        figure()
        plot(-dat.load(1:end-r)/(.25*pi*16.25e-3^2), hStore.^-1, 'r')
        xlabel("Pressure [Pa]")
        ylabel("Contact resistance h^{-1} [m^2 K W^{-1}]")
        ylim([0 1.5*max(prctile(hStore, [0 99.9]).^-1, [], "all")])
        grid on
        grid minor
    end
end

