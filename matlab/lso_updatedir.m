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
    % Form the dp / dgamma matrix.
    %

my_diag = @(x) spdiags(x(:), 0, numel(x), numel(x)); % Helper function.
dp_dg = repmat(my_diag(0.5 * sign(phi)), 1, 4);


    %
    % Form the dgamma / dphi matrix.
    %

[adj, on_border] = lso_priv_adjacents(phi); % Find cells on border.

% If no cells on border, return dphi = 0.
if isempty(find(on_border))
    dphi = 0 * dp;
    return
end

% This matrix connects how values of gamma are affected by values of phi.
dg_dphi = [my_dg_dphi(phi, adj{1}, 1); ...
    my_dg_dphi(phi, adj{2}, -1); ...
    my_dg_dphi(phi, adj{3}, dims(1)); ...
    my_dg_dphi(phi, adj{4}, -dims(1))];


    % 
    % Form selection matrix for the values of phi which are next to a border 
    % and active.
    %

ind = find(on_border & sel);
S_phi = sparse(1:length(ind), ind, ones(length(ind), 1), length(ind), N);


    %
    % Solve the following least-squares problem:
    %   minimize || A*x - b ||^2
    % 

A = dp_dg * dg_dphi * S_phi';
b = dp(:);

% When solving, add a small I since A is very often singular.
x = (A'*A + 1e-8 * speye(size(A,2))) \ (A' * b); 

% Insert values of x.
dphi = reshape(S_phi' * x, size(phi));



    %
    % Help construct the dg_dphi matrix.
    %

function [A] = my_dg_dphi(phi, adj, s)

N = numel(phi); % Number of elements in the grid.

% Find the cells whose gamma points are viable.
ind = find(adj);

% Form the submatrix.
A = sparse(ind, ind, -phi(ind+s)./(phi(ind)-phi(ind+s)).^2, N, N) + ...
    sparse(ind, ind+s, phi(ind)./(phi(ind)-phi(ind+s)).^2, N, N);
