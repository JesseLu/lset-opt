function [phi, err_bnd, err_p] = lso_regularize(phi_hat)
% [PHI, ERR_BND, ERR_P] = LSO_REGULARIZE(PHI_HAT)
%
% Description
%     Given proposed level-set function PHI_HAT, produces a "regularized"
%     level-set function PHI. The boundary points of PHI_HAT and PHI are 
%     identical. See the README or the theory manuscript for what is meant by
%     regularized.
% 
% Inputs
%     PHI_HAT: 2d array (level-set function).
%         The proposed level-set function.
% 
% Outputs
%     PHI: 2d array (level-set function).
%         The "regularized" level-set function.
% 
%     ERR_BND: Positive number.
%         The maximum displacement in either the x- or y-direction of a 
%         boundary point.
% 
%     ERR_P: Positive number.
%         The maximum change in the fractional filling.
% 
% Examples
%     % Generate random level-set and regularize.
%     phi_hat = rand(40, 30) - 0.2;
%     [phi, err_b, err_p] = lso_regularize(phi_hat);
%     lso_plot(phi);        


dims = size(phi_hat); % Size of the 2D grid.
N = prod(dims); % Number of elements in the grid.


    % 
    % Step 1: Set all phi_hat == 0 to phi_hat = eps.
    %

phi_hat = phi_hat + eps * (phi_hat == 0);


    %
    % Step 2: Fix cells which are not adjacent to boundaries.
    %

[adj, on_border] = lso_priv_adjacents(phi_hat); % Find border cells.
phi = (~on_border) .* sign(phi_hat); % Fix values of non-border cells.

% If there are no borders, we're done.
if isempty(find(on_border))
    err_bnd = 0;
    err_p = 0;
    return
end


    %
    % Step 3: Initialize remaining cells by solving minimization problem.
    %

A_ind = reshape(1:N, dims); % Reference matrix.

% Form difference matrix.
D = @(s) sparse(repmat(1:prod(dims-s), 1, 2), ... % For arbitrary differences.
    [A_ind(1:end-s(1),1:end-s(2)), A_ind(1+s(1):end,1+s(2):end)], ...
    [ones(dims-s), -ones(dims-s)], prod(dims-s), prod(dims));
D = [D([1 0]); D([0 1])]; % Concatenate both horizontal and vertical.

% Form selection matrix for cells on the border.
my_sel = @(ind) sparse(1:length(ind), ind, ones(length(ind),1), length(ind), N);
S_on = my_sel(find(on_border));

% Form equality constrained matrix.
my_eq = @(ind, shift) ... % Forms half of the matrix.
    sparse(repmat(1:length(ind), 1, 2), [ind, ind+shift], ...
    [phi_hat(ind+shift), -phi_hat(ind)], length(ind), N);
A = [my_eq(find(adj{1}), 1); my_eq(find(adj{3}), dims(1))];

% Solve for phi.
C = A * S_on';
b = sign(S_on * phi_hat(:));
d = zeros(size(A,1), 1);
x = lso_priv_solver(b, C, d);

% Form regularized phi.
phi = reshape(phi(:) + S_on'*x, dims);


    %
    % Step 4: Check answer (simple tests).
    %

% % Plot phi.
% subplot 121; lso_plot(phi);
% subplot 122; lso_plot(phi_hat);
% pause

% Make sure none of the values of phi has changed sign.
if any(sign(phi) ~= sign(phi_hat))
    error('During initialization, PHI changed sign!');
end

% Make sure boundary displacement error is acceptable.
[x_phi, y_phi] = lso_boundaries(phi);
[x_phi_hat, y_phi_hat] = lso_boundaries(phi_hat);
err_bnd = max([abs(x_phi - x_phi_hat); abs(y_phi - y_phi_hat)]);

if (err_bnd >= 1e-10)
    warning('Boundary displacement error, %e, exceeds threshold.', err_bnd);
end

% Make sure that the fractional-filling has not changed significantly.
err_p = max(max(abs(lso_fracfill(phi) - lso_fracfill(phi_hat))));
if (err_p >= 1e-10)
    warning('Fractional-filling error, %e, exceeds threshold.', err_p);
end

