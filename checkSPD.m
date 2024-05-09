function isSPD = checkSPD(A, tol)
    % checkSymmetry Checks if a matrix A is SPD within a specified tolerance.
    %   isSPD = checkSPD(A, tol) returns true if matrix A is SPD
    %   within the tolerance specified by tol. If tol is not specified, the
    %   function uses a default tolerance based on the machine precision.
    %
    %   Inputs:
    %       A - Matrix to be checked (can be dense or sparse)
    %       tol - Tolerance for comparing symmetry (optional)
    %
    %   Output:
    %       isSPD - Boolean flag indicating whether A is symmetric

    % Set default tolerance if not specified
    if nargin < 2 || isempty(tol)
      tol = max(size(A)) * eps(norm(A, 'fro'));
    end

    % can only checkPD since checkSymmetry is called in checkPD
    isSPD = checkPD(A, tol);
end

function isPD = checkPD(A, tol)
    % only symmetric matrix can be PD. 
    if ~checkSymmetry(A, tol)
        isPD = false;
        return;
    end

    % Attempt to use Cholesky decomposition to test positive definiteness
    try
        % Cholesky decomposition, which fails if A is not positive definite
        chol(A);
        isPD = true; % Cholesky successful, A is positive definite
    catch
        isPD = false; % Cholesky failed, A is not positive definite
    end
end

function isSymmetric = checkSymmetry(A, tol)
    % Check if A is equal to its transpose within the tolerance
    diffMatrix = abs(A - A');
    isSymmetric = all(diffMatrix(:) <= tol);

    %% can be removed (check different pairs)
    % Find pairs that are not exactly equal but within tolerance
    % [i, j] = find((diffMatrix > 0) & (diffMatrix <= tol));
    % diffPairs = zeros(length(i), 2);
    % 
    % for k = 1:length(i)
    %     if i(k) > j(k)  % Store only one of each pair since matrix is symmetric
    %         diffPairs(k, :) = [A(i(k), j(k)), A(j(k), i(k))];
    %     end
    % end
    % 
    % % Remove zero rows caused by the upper triangular check (i > j)
    % diffPairs( ~any(diffPairs, 2), : ) = [];
end