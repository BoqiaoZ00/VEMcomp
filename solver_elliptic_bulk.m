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
    RHS = M(:, bulknodes) * f(P(bulknodes, :)); 

    ubulk = LHS\RHS;
    u = zeros(length(M),1);
    u(bulknodes,:) = ubulk;

end