function [dphi] = lso_updatedir(phi, dp, sel)
% [DPHI] = LSO_UPDATEDIR(PHI, DP, SEL)
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
%     SEL: 2d array.
%         Marks the active (changeable) cells with a 1 and inactive cells with
%         0. 
%             
% Outputs
%     DPHI: 2d array.
%         The derivative in the level-set function.


dims = size(phi); % Size of grid.
N = prod(dims); % Number of elements in the grid.


    % 
    % Form selection matrix for the values of phi which are next to a border 
    % and active.
    %

[adj, on_border] = lso_priv_adjacents(phi); % Find cells on border.
ind = find(on_border & sel); % Indices of active, on-border cells.

% If no active on-border cells, return dphi = 0.
if isempty(ind)
    dphi = 0 * dp;
    return
end
 
% Form sparse selection matrix, which "picks out" the relevant cells.
S_phi = sparse(1:length(ind), ind, ones(length(ind), 1), length(ind), N);


    %
    % Form the dp / dphi matrix.
    %

g = lso_priv_gamma(phi); % Find gamma values (positive distances to boundary).

% Determine index of the minimum-distance partner cell.
[gamma, dir] = min([g{1}(ind), g{2}(ind), g{3}(ind), g{4}(ind)], [], 2);
ind_i = ind + ((dir == 1) - (dir == 2)) + dims(1) * ((dir == 3) - (dir == 4));

% Form the matrix.
grad = abs(phi(ind) - phi(ind_i)); % Gradient.
% dp_dphi = sparse(ind, ind, abs(phi(ind_i)) ./ grad.^2, N, N) + ...
%     sparse(ind, ind_i, abs(phi(ind)) ./ grad.^2, N, N);
dp_dphi = sparse(ind, ind, abs(phi(ind_i)) ./ grad.^2, N, N);


    %
    % Solve the following least-squares problem:
    %   minimize || A*x - b ||^2
    % 

A = dp_dphi * S_phi';
b = dp(:);

% When solving, add a small I since A is very often singular.
% x = A \ b;
x = (A'*A + 1e-10 * speye(size(A,2))) \ (A' * b); 

% Insert values of x.
dphi = reshape(S_phi' * x, size(phi));


