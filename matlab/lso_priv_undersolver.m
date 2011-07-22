function [x] = lso_priv_undersolver(x0, A, b)
% X = LSO_PRIV_UNDERSOLVER(X0, A, B)
% 
% Description
%     Solves the following problem for variable x:
%         minimize ||x - x0||^2
%         subject to A * x == b
%     For the special case where A is fat and the rows of A are not 
%     linearly independent.
% 
%     This is accomplished by taking the QR decomposition of A^T, in order to
%     eliminate linearly dependent rows A. The following formula is then 
%     applied:
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
    % Verify inputs.
    %

if (~issparse(A)) % Make sure A is sparse.
    error('Matrix A must be sparse.');
end 

% if (size(A,1) > size(A,2)) % Make sure A is fat.
%     error('Matrix A must be fat (size of A is %d x %d).', size(A));
% end


    % 
    % Take QR-decomposition of A and eliminate linearly-dependent rows of A,
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
R = R(:,ind);
R = R(1:length(ind),:);
A = A(ind,:);
b = b(ind);



    %
    % Solve for x.
    %

y = R \ (R' \ (A * x0 - b));
x = x0 - A' * y;


    % 
    % Calculate error.
    %

