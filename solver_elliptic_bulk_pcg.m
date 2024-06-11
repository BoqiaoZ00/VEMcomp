function u = solver_elliptic_bulk_pcg(D,alpha,f,P,M,K,R,bcond)
    
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

   
    maxit = 500;
    tol = 1e-3;
    pre = ichol(LHS);
    LHS = 0.5 * (LHS + LHS');
    ubulk = pcg(LHS, RHS, tol, maxit, pre, pre');

    u = zeros(length(M),1);
    u(bulknodes,:) = ubulk;

end
