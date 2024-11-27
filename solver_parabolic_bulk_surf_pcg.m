function [u, v, t, uprime_norm, vprime_norm, u_average, v_average] = solver_parabolic_bulk_surf_pcg(dOmega, dGamma, f, g, h, P, M, MS, K, KS, R, T, tau, u0, v0)

n = length(dOmega);
m = length(dGamma);
PGamma = R' * P;

% Precompute LHS matrices and preconditioners
LHS_Omega = cell(1, n);
M_precond_Omega = cell(1, n);
for i = 1:n
    LHS_Omega{i} = M + tau * dOmega(i) * K; 
    M_precond_Omega{i} = ichol(LHS_Omega{i}, struct('type', 'ict', 'droptol', 1e-3));
end
LHS_Gamma = cell(1, m);
M_precond_Gamma = cell(1, m);
for i = 1:m
    LHS_Gamma{i} = MS + tau * dGamma(i) * KS; 
    M_precond_Gamma{i} = ichol(LHS_Gamma{i}, struct('type', 'ict', 'droptol', 1e-3)); 
end

NT = ceil(T / tau);
u = u0;
v = v0;
if nargout >= 4
    % spatial L2 norm of time derivative of fist component of u and v
    uprime_norm = zeros(1, NT);
    vprime_norm = zeros(1, NT);
end
if nargout >= 6
    % spatial average of first component of u and v
    u_average = zeros(1, NT);
    v_average = zeros(1, NT);
end
progress_handle = waitbar(0);
axes_handle = findobj(progress_handle, 'type', 'axes');
title_handle = get(axes_handle, 'title');
set(title_handle, 'FontSize', 18);
waitbar(0, progress_handle, 'Timestepping in progress: 0 %')
percent_prev = 0;

tol = 1e-10;
maxit = 1000;

for i = 0:NT-1

    percent_new = round(i * 100 / NT);
    if percent_new > percent_prev
        waitbar(i / NT, progress_handle, sprintf('Timestepping in progress: %d%%', percent_new));
        percent_prev = percent_new;
    end

    RHS_Omega = zeros(size(M, 1), n);
    for j = 1:n
        RHS_Omega(:, j) = M * (u(:, j) + tau * f{j}(u, P, i * tau)) + tau * R * MS * h{j}(R' * u, v, PGamma, i * tau); 
    end
    RHS_Gamma = zeros(size(MS, 1), m);
    for j = 1:m
        RHS_Gamma(:, j) = MS * (v(:, j) + tau * g{j}(R' * u, v, PGamma, i * tau)); 
    end

    unew = zeros(size(LHS_Omega{1}, 1), n);
    for j = 1:n
        [unew(:, j), flag] = pcg(LHS_Omega{j}, RHS_Omega(:, j), tol, maxit, M_precond_Omega{j}, M_precond_Omega{j}', u(:, j)); 
        if flag ~= 0
            warning('PCG did not converge for component %d at time step %d', j, i);
        end
    end
    vnew = zeros(size(LHS_Gamma{1}, 1), m);
    for j = 1:m
        [vnew(:, j), flag] = pcg(LHS_Gamma{j}, RHS_Gamma(:, j), tol, maxit, M_precond_Gamma{j}, M_precond_Gamma{j}', v(:, j)); 
        if flag ~= 0
            warning('PCG did not converge for component %d at time step %d', j, i);
        end
    end

    if nargout >= 4
        incr_u = unew(:, 1) - u(:, 1);
        incr_v = vnew(:, 1) - v(:, 1);
        uprime_norm(i + 1) = incr_u' * M * incr_u;
        vprime_norm(i + 1) = incr_v' * MS * incr_v;
    end
    if nargout >= 6
        u_average(i + 1) = sum(M * unew(:, 1));
        v_average(i + 1) = sum(MS * vnew(:, 1));
    end

    u = unew;
    v = vnew;
end

close(progress_handle);

if nargout >= 3
    t = linspace(tau, NT * tau, NT);
end
if nargout >= 4
    uprime_norm = sqrt(uprime_norm);
    vprime_norm = sqrt(vprime_norm);
end
if nargout >= 6
    u_average = u_average / sum(sum(M));
    v_average = v_average / sum(sum(MS));
end

end
