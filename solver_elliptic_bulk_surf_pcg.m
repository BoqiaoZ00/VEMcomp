function  [u, v]  = solver_elliptic_bulk_surf_pcg(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, varargin)
    % Check the number of input arguments
    numArgs = nargin;
    
    switch numArgs
         case 1
            tol = varargin{1};
            [u, v] = solver_elliptic_bulk_surf_pcg_tol(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol);
        case 2
            tol = varargin{1};
            maxit = varargin{2};
             [u, v] = solver_elliptic_bulk_surf_pcg_tolmaxit(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol, maxit);
         case 3
            tol = varargin{1};
            maxit = varargin{2};
            pcg_M =  varargin{3};
            [u, v] = solver_elliptic_bulk_surf_pcg_tolmaxitM(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol, maxit, pcg_M);
        case 4
            tol = varargin{1}; 
            maxit = varargin{2};
            M1 = varargin{3};
            M2 = varargin{4};
            [u, v]  = solver_elliptic_bulk_surf_pcg_M1M2(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol, maxit, M1, M2);
         case 5
            tol = varargin{1}; 
            maxit = varargin{2};
            M1 = varargin{3};
            M2 = varargin{4};
            x0 = varargin{5};
            [u, v] = solver_elliptic_bulk_surf_pcg_tolmaxitM1M2x0(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol, maxit, M1, M2, x0);
        otherwise
            error('Invalid number of input arguments');
    end
end

function [u, v] = solver_elliptic_bulk_surf_pcg_tol(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol)

N = length(M);
LHS = [DOmega * K + alpha * M - gamma * (R * MS * R'), -delta * (R * MS);
       sparse(size(MS, 1), N), DGamma * KS + beta * MS];

RHS = [M * f(P); MS * g(R' * P)];

% Solve using PCG
x = pcg(LHS, RHS, tol);

% Extract solutions
u = x(1:N, 1);
v = x(N + 1:end, 1);

end

function [u, v] = solver_elliptic_bulk_surf_pcg_tolmaxit(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol, maxit)

N = length(M);
LHS = [DOmega * K + alpha * M - gamma * (R * MS * R'), -delta * (R * MS);
       sparse(size(MS, 1), N), DGamma * KS + beta * MS];

RHS = [M * f(P); MS * g(R' * P)];

% Solve using PCG
x = pcg(LHS, RHS, tol, maxit);

% Extract solutions
u = x(1:N, 1);
v = x(N + 1:end, 1);

end

function [u, v] = solver_elliptic_bulk_surf_pcg_tolmaxitM(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol, maxit, pcg_M)

N = length(M);
LHS = [DOmega * K + alpha * M - gamma * (R * MS * R'), -delta * (R * MS);
       sparse(size(MS, 1), N), DGamma * KS + beta * MS];

RHS = [M * f(P); MS * g(R' * P)];

% Solve using PCG
x = pcg(LHS, RHS, tol, maxit, pcg_M);

% Extract solutions
u = x(1:N, 1);
v = x(N + 1:end, 1);

end

function [u, v] = solver_elliptic_bulk_surf_pcg_tolmaxitM1M2(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol, maxit, M1, M2)

N = length(M);
LHS = [DOmega * K + alpha * M - gamma * (R * MS * R'), -delta * (R * MS);
       sparse(size(MS, 1), N), DGamma * KS + beta * MS];

RHS = [M * f(P); MS * g(R' * P)];

% Use the default initial guess of PCG
x0 = zeros(size(RHS));

% Solve using PCG
x = pcg(LHS, RHS, tol, maxit, M1, M2, x0);

% Extract solutions
u = x(1:N, 1);
v = x(N + 1:end, 1);

end

function [u, v] = solver_elliptic_bulk_surf_pcg_tolmaxitM1M2x0(DOmega, DGamma, alpha, beta, f, g, P, M, MS, K, KS, R, gamma, delta, tol, maxit, M1, M2, x0)

N = length(M);
LHS = [DOmega * K + alpha * M - gamma * (R * MS * R'), -delta * (R * MS);
       sparse(size(MS, 1), N), DGamma * KS + beta * MS];

RHS = [M * f(P); MS * g(R' * P)];

% Solve using PCG
x = pcg(LHS, RHS, tol, maxit, M1, M2, x0);

% Extract solutions
u = x(1:N, 1);
v = x(N + 1:end, 1);

end