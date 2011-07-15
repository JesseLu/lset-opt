function [x, v, solve_time] = lso_priv_quadeq (A, b, C, d)
% [X, V, SOLVE_TIME] = LSO_PRIV_QUADEQ (A, B, C, D)
% 
% Description
%     Solve the quadratic minimization problem with linear equality given by
%         minimize ||Ax - b||^2
%         subject to C x = d.

n = size(A, 1);
p = size(C, 1);

b = A'*b;
A = A'*A;

% Solve.
tic;
w = A \ [C', b];
v = (C * w(:,1:p)) \ (C * w(:,end) - d);
x = w(:,end) - w(:,1:p) *  v;
solve_time = toc;

% % Alternative: directly solve the compound matrix.
% tic
% z = [A, C; C', zeros(size(C,2))] \ [b; d];
% x = z(1:n);
% v = z(n+1:end);
% solve_time = toc;
