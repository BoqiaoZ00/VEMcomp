function [v, t, vprime_norm, v_average] = solver_parabolic_surf_pcg(D,g,P,MS,KS,R,T,tau,v0)
    
    n = length(D);
    PS = R'*P;
    NT = ceil(T/tau);

    v = v0;
    v_new = zeros(size(v0));
    if nargout >= 3
        vprime_norm = zeros(1,NT);
    end
    if nargout >= 4
        v_average = zeros(1,NT);
    end

    size(KS)

    LHS = cell(1, n);
    M = cell(1, n);
    for j=1:n
        LHS{j} = tau*D(j)*KS + MS; %#ok
        M{j} = ichol(LHS{j}, struct('type', 'ict', 'droptol', 1e-3)); %#ok
    end

    for i=0:NT-1
        for j=1:n
            RHS = MS*(v(:,j) + tau*g{j}(v, PS, i*tau));


            % Use PCG to solve the linear system
            % Assuming LHS is SPD; if not, additional work may be needed to ensure it
            tol = 1e-6;
            maxit = 50;
            [v_new(:,j), flag] = pcg(LHS{j}, RHS, tol, maxit, M{j}, M{j}', v(:,j)); 
            if flag ~= 0
                error('PCG did not converge');
            end
        end
        if nargout >= 3
            incr = v_new(:,1) - v(:,1);
            vprime_norm(i+1) = incr'*MS*incr;
        end
        if nargout >= 4
            v_average(i+1) = sum(MS*v_new(:,1));
        end
        v = v_new;
    end

    if nargout >= 2
        t = linspace(tau, NT*tau, NT);
    end
    if nargout >= 3
        vprime_norm = sqrt(vprime_norm);
    end
    if nargout >= 4
        v_average = v_average / sum(sum(MS));
    end

end