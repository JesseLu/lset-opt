function [dphi] = lso_updatedir(phi, dp)
% [DPHI] = LSO_UPDATEDIR(PHI, DP)
% 
% Description
%     Find a derivative in phi to approximate the derivative in p.
% 
% Inputs
%     PHI: 2d array (level-set function).
% 
%     DP: 2d array.
%         The derivative in the fractional-filling (p).
%             
% Outputs
%     DPHI: 2d array.
%         The derivative in the level-set function.

dims = size(phi); % Size of grid.
N = prod(dims); % Number of elements in the grid.
 

    %
    % Form the dp / dgamma matrix.
    %

gamma = lso_priv_gamma(phi); % Find boundary points.
width = 2 * sign(phi) .* (gamma{1} + gamma{2});
height = 2 * sign(phi) .* (gamma{3} + gamma{4});

my_diag = @(x) spdiags(x(:), 0, numel(x), numel(x)); % Helper function.
dp_dg = [my_diag(height), my_diag(height), my_diag(width), my_diag(width)];


    %
    % Form the dgamma / dphi matrix.
    %

[adj, on_border] = lso_priv_adjacents(phi); % Find cells on border.

% This matrix connects how values of gamma are affected by values of phi.
dg_dphi = [my_dg_dphi(phi, adj{1}, 1); ...
    my_dg_dphi(phi, adj{2}, -1); ...
    my_dg_dphi(phi, adj{3}, dims(1)); ...
    my_dg_dphi(phi, adj{4}, -dims(1))];


    % 
    % Form selection matrix for the values of phi which are next to a border.
    %

ind = find(on_border);
S_phi = sparse(1:length(ind), ind, ones(length(ind), 1), length(ind), N);


    %
    % Solve the following problem:
    %     minimize ||dphi||^2
    %     subject to C dphi == dp
    % 

C = dp_dg * dg_dphi * S_phi';
d = dp(:);
b = zeros(size(C,2), 1);
x = b - C' * ((C * C') \ (C * b - d));
dphi = reshape(S_phi' * x, size(phi));



    %
    % Help construct the dg_dphi matrix.
    %

function [A] = my_dg_dphi(phi, adj, s)

N = numel(phi); % Number of elements in the grid.

% Find the cells whose gamma points are viable.
ind = find(adj);
if (s > 0) % Break ties.
    ind = ind(find(abs(phi(ind)) < abs(phi(ind+s))));
else
    ind = ind(find(abs(phi(ind)) <= abs(phi(ind+s))));
end

% Form the submatrix.
A = sparse(ind, ind, -phi(ind+s)./(phi(ind)-phi(ind+s)).^2, N, N) + ...
    sparse(ind, ind+s, phi(ind)./(phi(ind)-phi(ind+s)).^2, N, N);
