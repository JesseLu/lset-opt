function [x, v, solve_time] = la_quadeq (A, b, C, d)

n = size(A, 1);
p = size(C, 2);

% Solve.
tic;
w = A \ [C, b];
v = (C' * w(:,1:p)) \ (C' * w(:,end) - d);
x = w(:,end) - w(:,1:p) *  v;
solve_time = toc;

% % Alternative: directly solve the compound matrix.
% tic
% z = [A, C; C', zeros(size(C,2))] \ [b; d];
% x = z(1:n);
% v = z(n+1:end);
% solve_time = toc;
