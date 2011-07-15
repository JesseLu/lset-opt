function [phi] = lso_initialize(phi_hat)

dims = size(phi_hat); % Size of the 2D grid.
N = prod(dims);


    % 
    % Step 1: Set all phi_hat == 0 to phi_hat = eps.
    %

phi_hat = phi_hat + eps * (phi_hat == 0);


    %
    % Step 2: Find cells which are adjacent to a boundary point.
    %

% Boundary present to the right of the cell
adj_right = sign(phi_hat(1:end-1,:)) ~= sign(phi_hat(2:end,:));

% Boundary present downward of the cell
adj_down = sign(phi_hat(:,1:end-1)) ~= sign(phi_hat(:,2:end));

% Find which cells are next to a boundary point.
on_border = [adj_right; zeros(1, dims(2))] | ...
    [zeros(1, dims(2)); adj_right] | ...
    [adj_down, zeros(dims(1), 1)] | ...
    [zeros(dims(1), 1), adj_down];

% Make the adj_* arrays the same dimensions as phi_hat.
adj_right = cat(1, adj_right, zeros(1, dims(2)));
adj_down = cat(2, adj_down, zeros(dims(1), 1));


    %
    % Step 2: Fix cells which are not adjacent to boundaries.
    %

phi = (~on_border) .* sign(phi_hat);


    %
    % Step 3: Initialize remaining cells by solving minimization problem.
    %

A_ind = reshape(1:N, dims); % Reference matrix.

% Form difference matrix.
D = @(s) sparse(repmat(1:prod(dims-s), 1, 2), ... % For arbitrary differences.
    [A_ind(1:end-s(1),1:end-s(2)), A_ind(1+s(1):end,1+s(2):end)], ...
    [ones(dims-s), -ones(dims-s)], prod(dims-s), prod(dims));
D = [D([1 0]); D([0 1])]; % Concatenate both horizontal and vertical.

% Form selection matrices for cells on and off the border.
my_sel = @(ind) sparse(1:length(ind), ind, ones(length(ind),1), length(ind), N);
S_on = my_sel(find(on_border));
S_off = my_sel(find(~on_border));

% Form equality constrained matrix.
my_eq = @(ind, shift) ... % Forms half of the matrix.
    sparse(repmat(1:length(ind), 1, 2), [ind; ind+shift], ...
    [phi_hat(ind+shift); -phi_hat(ind)], length(ind), N);
A = [my_eq(find(adj_right), 1); my_eq(find(adj_down), dims(1))];

% Solve for phi.
% spy(S_on*A')
% A*S_on'*S_on*phi_hat(:)
% A*phi_hat(:)
% S_on
% imagesc(reshape(S_on'*S_on*phi_hat(:), dims)');
% hold on; spy(adj_right);
% phi_hat
% adj_right

[x, v, solve_time] = la_quadeq(S_on*D'*D*S_on', -S_on*D'*D*phi(:), ...
    S_on*A', zeros(size(A,1), 1));

phi = reshape(phi(:) + S_on'*x, dims);



% Plot phi.
imagesc(phi');
axis equal tight;
phi
phi_hat
