function [u, v] = solver_elliptic_bulk_surf_pcg(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta)

N = length(M);
LHS = [DOmega * K + alpha * M - gamma * (R * MS * R'), -delta * (R * MS);
       sparse(size(MS, 1), N), DGamma * KS + beta * MS];

RHS = [M * f(P); MS * g(R' * P)];

% Preconditioner using Incomplete Cholesky factorization
try
    M_precond = ichol(LHS, struct('type', 'ict', 'droptol', 1e-3));
catch
    warning('ichol failed with droptol=1e-3, trying with droptol=1e-2');
    try
        M_precond = ichol(LHS, struct('type', 'ict', 'droptol', 1e-2));
    catch
        warning('ichol failed with droptol=1e-2, trying with droptol=1e-1');
        M_precond = ichol(LHS, struct('type', 'ict', 'droptol', 1e-1));
    end
end

% Initial guess
x0 = diag(diag(ones(size(RHS))));

size(LHS)

% Tolerance and maximum iterations for PCG
tol = 1e-4;
maxit = 3000;

% Solve using PCG
x = pcg(LHS, RHS, tol, maxit, M_precond, M_precond', x0);

% Extract solutions
u = x(1:N, 1);
v = x(N + 1:end, 1);

end