function [x, y] = lso_boundaries(phi)
% [X, Y] = LSO_BOUNDARIES(PHI)
% 
% Description
%     Compute a list of boundary points (X, Y) from the level-set function PHI.
% 
% Inputs
%     PHI: 2d array (level-set function).
% 
% Outputs
%     X, Y: Vectors.
%         The elements of X and Y denote the location of the boundary points
%         implicitly described by PHI. The order of the points is unique in 
%         the sense that any level-set functions describing the same boundary 
%         points will produce identically ordered (X, Y).

dims = size(phi); % Size of grid.

% Horizontal and vertical borders.
adj_right = cat(1, sign(phi(1:end-1,:)) ~= sign(phi(2:end,:)), ...
    zeros(1, dims(2)));
adj_down = cat(2, sign(phi(:,1:end-1)) ~= sign(phi(:,2:end)), ...
    zeros(dims(1), 1));

% Calculate horizontal zero-crossing points.
ind = find(adj_right);
[hx, hy] = ind2sub(dims, ind);
hx = hx + phi(ind) ./ (phi(ind) - phi(ind+1));

% Calculate vertical zero-crossing points.
ind = find(adj_down);
[vx, vy] = ind2sub(dims, ind);
vy = vy + phi(ind) ./ (phi(ind) - phi(ind+dims(1)));

% Concatenate lists.
x = [hx(:); vx(:)];
y = [hy(:); vy(:)];

% % Plot.
% imagesc(phi'); hold on; plot(x, y, 'r.'); axis equal tight;
