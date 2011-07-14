function [phi] = lso_initialize(phi_hat)

dims = size(phi_hat); % Size of the 2D grid.


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


    %
    % Step 2: Fix cells which are not adjacent to boundaries.
    %

phi = (~on_border) .* sign(phi_hat);


    %
    % Step 3: Initialize remaining cells by solving minimization problem.
    %

% Form difference matrix.

% Form selection matrices for cells on and off the border.

% Form equality constrained matrix.

% Solve for phi.


% Plot phi.
imagesc(phi');
axis equal tight;
