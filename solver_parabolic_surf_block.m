function [v, t, vprime_norm, v_average] = solver_parabolic_surf_block(D,g,P,MS,KS,R,T,tau,v0)
    
    n = length(D);
    PS = R'*P;
    NT = ceil(T/tau);

    v = v0;
    v_block = reshape(v, [], 1);
    if nargout >= 3
        vprime_norm = zeros(1,NT);
    end
    if nargout >= 4
        v_average = zeros(1,NT);
    end

    LHS = cell(1, n);
    for j=1:n
        LHS{j} = tau*D(j)*KS + MS; 
    end

    LHS_block = blkdiag(LHS{:});
    %precond = diag(diag(LHS_block));
    precond = ichol(LHS_block, struct('type', 'ict', 'droptol', 1e-3)); 
    tol = 1e-6;
    maxit = 200;
    p = size(MS, 1);

    MS_cell = repmat({MS}, 1, n);
    MS_block = blkdiag(MS_cell{:});
    % [r, ~] = size(v);
    % g_block = zeros(n*r, 1);

    tic
    for i=0:NT-1
        % index = 1;
        % for j = 1:n
        %     g_block(index:index+r-1) = g{j}(v, PS, i * tau);
        %     index = index + r;
        % end
        RHS_block = MS_block * v_block; % + tau * g_block);
        [v_new_block, flag] = pcg(LHS_block, RHS_block, tol, maxit, precond, precond', v_block); 
        if flag ~= 0
            error('PCG did not converge');
        end
        % index=1;
        % for j = 1:n
        %     v(:, j) = v_block(index:index + p - 1);
        %     index = index + p;
        % end
        % if nargout >= 3
        %     v_new = reshape(v_new_block, [], n);
        %     incr = v_new(:,1) - v(:,1);
        %     vprime_norm(i+1) = incr'*MS*incr;
        % end
        % if nargout >= 4
        %     v_new = reshape(v_new_block, [], n);
        %     v_average(i+1) = sum(MS*v_new(:,1));
        % end
        v_block = v_new_block;
    end
    toc

    index=1;
    for j = 1:n
        v(:, j) = v_block(index:index + p - 1);
        index = index + p;
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