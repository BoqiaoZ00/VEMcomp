function [v] = solver_elliptic_surf_pcg(D, alpha, g, P, MS, KS, R)

    LHS = D*KS + alpha*MS;
    RHS = MS*g(R'*P);

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
    
    % Tolerance and maximum iterations for PCG
    tol = 1e-4;
    maxit = 3000;
    
    % Solve using PCG
    v = pcg(LHS, RHS, tol, maxit, M_precond, M_precond', x0);
    
    % Check for convergence
    % if flag ~= 0
    %     error('PCG did not converge');
    % end
    

end