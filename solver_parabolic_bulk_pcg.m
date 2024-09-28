function [u, t, uprime_norm, u_average] = solver_parabolic_bulk_pcg(D, f, P, M, K, R, bcond, T, tau, u0)
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
    M_precond = cell(1, n);
    for i = 1:n
        LHS{i} = tau * D(i) * Kbcond + Mbcond;
        M_precond{i} = ichol(LHS{i}, struct('type', 'ict', 'droptol', 1e-3));
    end

    NT = ceil(T / tau);
    u = u0;
    ubcond = zeros(length(bulknodes), n);  % Preallocate ubcond
    RHS = zeros(length(bulknodes), n);  % Preallocate RHS

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

        for j = 1:n
            % Compute the RHS for bulk nodes only
            RHS(:, j) = M(bulknodes, :) * (u(:, j) + tau * f{j}(u, P, i * tau));
        end

        for j = 1:n
            tol = 1e-10;
            maxit = 100;
            % Use u(bulknodes, j) as the initial guess for the PCG solver
            [ubcond(:, j), flag] = pcg(LHS{j}, RHS(:, j), tol, maxit, M_precond{j}, M_precond{j}', u(bulknodes, j));
            if flag ~= 0
                warning('PCG did not converge for component %d at step %d', j, i + 1);
            end
        end

        u_new = zeros(length(M), n);  % Initialize u_new here
        u_new(bulknodes, :) = ubcond;

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