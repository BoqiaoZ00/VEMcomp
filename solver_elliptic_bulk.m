function u = solver_elliptic_bulk(D,alpha,f,P,M,K,R,bcond)
    
    bulknodes = (1:length(M))';
    if strcmp(bcond, 'dir')
        boundarynodes = sum(R,2) == 1;
    else
        boundarynodes = [];
    end
    bulknodes(boundarynodes) = [];

    Mdir = M(bulknodes, bulknodes);
    Kdir = K(bulknodes, bulknodes);
    
    LHS = D*Kdir + alpha*Mdir;
    rhs = M*f(P);
    RHS = rhs(bulknodes);

    % start
    tic 
    ubulk = LHS\RHS;
    toc
    

    tic
    % if (~checkSPD(LHS, 1e-12))
    %     error("cannot use PCG")
    % end
    % TODO: prove mathematically that LHS is SPD
    if (~isequal(LHS, LHS'))
        % LHS is almost (but numerically not) symmetric positive definite
        LHS = 0.5 * (LHS + LHS');
    end
    maxit = 500;
    tol = 1e-10;
    pre = ichol(LHS);
    utest = pcg(LHS, RHS, tol, maxit, pre, pre');
    toc
    % end
    
    u = zeros(length(M),1);
    u(bulknodes,:) = ubulk;

end