function [x, err] = lso_priv_solver(x0, A, b)
% [X, ERR] = LSO_PRIV_SOLVER(X0, A, B)
% 
% Description
%     Solves the following problem for variable x:
%         minimize ||x - x0||^2
%         subject to A * x == b
%     For the special case where the rows of A are not linearly independent.
% 
%     This is accomplished by taking the QR decomposition of A^T, in order to
%     eliminate linearly dependent rows A. The following formula is then 
%     applied (Ref. Boyd EE 263 Notes (Stanford) pg. 8-15):
%         x = x0 - A^T (A A^T)^-1 (A x0 - b),
%     where A and b have been modified to eliminate linearly dependent rows and
%     their corresponding entries.
% 
% Inputs
%     X0: Column vector (n x 1).
% 
%     A: Sparse matrix (m x n).
% 
%     B: Column vector (m x 1).
% 
% Outputs
%     X: Column vector (n x 1).
%         X is guaranteed to satisfy Ax = b (within numerical error) while 
%         minimizing ||x - x0||.
% 
%     ERR: Positive number.
%         The maximum element of |A*x - b|. If ERR > 1e-10, then a warning is
%         issued.

    %
    % Verify inputs.
    %

if (~issparse(A)) % Make sure A is sparse.
    error('Matrix A must be sparse.');
end 


    % 
    % Take QR-decomposition of A^T and eliminate linearly-dependent rows of A,
    % as well as corresponding entries in b.
    %

R = qr(A', 0); % Do not compute Q, since it is not needed.

% Find linearly-independent columns of R.
[i, j] = ind2sub(size(R), find(R));
ind = [];
for k = 1 : max(i)
    ind(end+1) = min(j(find(i == k)));
end

% Keep only linearly-independent entries.
R = R(:,ind); % Eliminate linearly dependent columns of R.
R = R(1:length(ind),:); % Remove rows in R.
A = A(ind,:); % Eliminate dependent rows of A.
b = b(ind); % Eliminate corresponding entries of b.


    %
    % Solve for x.
    %

% Here, the inverse of AA^T is replaced by R^T*R, since we have R already.
% Note: The use of R's may induce additional numerical errors. If these
% errors are unacceptable, there should be additional (cheap) computations
% to obtain additional accuracy in x.
x = x0 - A' * (R \ (R' \ (A * x0 - b)));


    % 
    % Calculate error.
    %

err = max(abs(A * x - b));
if (err > 1e-10)
    warning('Error in solver exceeds threshold (%e).', err);
end
