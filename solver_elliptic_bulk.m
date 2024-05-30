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
    % custom_spy(LHS);
    % spy(LHS);
    tic 
    ubulk = LHS\RHS;
    toc
    

    tic
    % % if (~checkSPD(LHS, 1e-12))
    % %     error("cannot use PCG")
    % % end
    % % TODO: prove mathematically that LHS is SPD
    % if (~isequal(LHS, LHS'))
    %     % LHS is almost (but numerically not) symmetric positive definite
    %     LHS = 0.5 * (LHS + LHS');
    % end
    maxit = 500;
    tol = 1e-6; % default value
    pre = ichol(LHS);
    LHS = 0.5 * (LHS + LHS');
    utest = pcg(LHS, RHS, tol, maxit, pre, pre');
    toc

    tic
    maxit = size(LHS,1);
    pre = ichol(LHS);
    LHS = 0.5 * (LHS + LHS');
    utest2 = gmres(LHS, RHS, [], tol, maxit, pre, pre');
    toc

    size(LHS)

    u = zeros(length(M),1);
    u(bulknodes,:) = ubulk;

end

% helper only used for testing
function custom_spy(A)
    [i, j, s] = find(A);  % Extract row indices, column indices, and values of non-zero entries
    figure;
    scatter(j, i, 20, abs(s), 'filled');  % Create a scatter plot with size and color based on the magnitude
    colormap('jet');  % Use a colormap to indicate magnitude
    colorbar;  % Show a colorbar as a legend for the magnitudes
    set(gca, 'YDir', 'reverse');  % Reverse y-axis to match the 'spy' view
    axis equal tight;  % Set axes to equal and tight
    xlabel('Column index');
    ylabel('Row index');
    title('Custom Spy Plot with Magnitude Information');
end