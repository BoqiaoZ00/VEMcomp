function [u, t, uprime_norm, u_average] = solver_parabolic_bulk_block(D, f, P, M, K, R, bcond, T, tau, u0)
    n = length(D);

    bulknodes = (1:length(M))';
    if strcmp(bcond, 'dir')
        boundarynodes = sum(R, 2) == 1;
    else
        boundarynodes = [];
    end
    bulknodes(boundarynodes) = [];

    Mbcond = sparse(M(bulknodes, bulknodes));
    Kbcond = sparse(K(bulknodes, bulknodes));

    % Precompute LHS matrices and preconditioners
    LHS = cell(1, n);
    for i = 1:n
        LHS{i} = tau * D(i) * Kbcond + Mbcond;
    end
    LHS_block = blkdiag(LHS{:});
    M_precond = ichol(LHS_block, struct('type', 'ict', 'droptol', 1e-3)); 

    NT = ceil(T / tau);
    u = u0;
    ubcond = zeros(length(bulknodes), n);  % Preallocate ubcond
    ubcond_block = zeros(length(bulknodes) * n, 1);  % Preallocate ubcond_block
    ubcond_block(:)=u(:);
    RHS = zeros(length(bulknodes), n);  % Preallocate RHS
    RHS_block = zeros(length(bulknodes) * n, 1);  % Preallocate RHS_block

    u_new = zeros(length(M), n); 
    u_initial_guess = zeros(size(bulknodes, 1) * n, 1);

    tol = 1e-10;
    maxit = 100;

    if nargout >= 3
        uprime_norm = zeros(1, NT);
    end
    if nargout >= 4
        u_average = zeros(1, NT);
    end

    progress_handle = waitbar(0);
    axes_handle = findobj(progress_handle, 'type', 'axes');
    title_handle = get(axes_handle, 'title');
    set(title_handle, 'FontSize', 18);
    waitbar(0, progress_handle, 'Timestepping in progress: 0 %')

    percent_prev = 0;
    for i = 0:NT-1
        percent_new = round(i * 100 / NT);
        if percent_new > percent_prev
            waitbar(i / NT, progress_handle, sprintf('Timestepping in progress: %d%%', percent_new));
            percent_prev = percent_new;
        end

        u(:)=ubcond_block(:);
        for j = 1:n
            % Compute the RHS for bulk nodes only
            RHS(:, j) = M(bulknodes, :) * (u(:, j) + tau * f{j}(u, P, i * tau));
        end
        RHS_block(:) = RHS(:);

        u_initial_guess(:) = u(bulknodes, :);

        [ubcond_block, flag] = pcg(LHS_block, RHS_block, tol, maxit, M_precond, M_precond', u_initial_guess);
        if flag ~= 0
            warning('PCG did not converge at step %d', i + 1);
        end

        u_new(bulknodes, :) = reshape(ubcond_block, size(bulknodes,1), n);

        if nargout >= 3
            incr = u_new(:, 1) - u(:, 1);
            uprime_norm(i + 1) = sqrt(incr' * M * incr);
        end
        if nargout >= 4
            u_average(i + 1) = sum(M * u_new(:, 1)) / sum(sum(M));
        end
        u = u_new;
    end

    close(progress_handle);

    if nargout >= 2
        t = linspace(tau, NT * tau, NT);
    end
    if nargout >= 3
        uprime_norm = sqrt(uprime_norm);
    end
    if nargout >= 4
        u_average = u_average / sum(sum(M));
    end
end