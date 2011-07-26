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

my_diag = @(x) spdiags(x(:), 0, numel(x), numel(x)); % Helper function.
dp_dg = repmat(my_diag(sign(phi)), 1, 4);


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
    % Form selection matrix for the values of phi which are next to a border.
    %

ind = find(on_border);
S_phi = sparse(1:length(ind), ind, ones(length(ind), 1), length(ind), N);


    %
    % Find weights for phi on border, but without control of boundary points.
    %

gamma = lso_priv_gamma(phi);
s = lso_priv_shifted(phi);
w_ind = ind(find(((phi(ind) .* dp(ind)) < 0) .* ...
    (gamma{1}(ind) == 0.5) .* (gamma{2}(ind) == 0.5) .* ...
    (gamma{3}(ind) == 0.5) .* (gamma{4}(ind) == 0.5))); 
phi_i = sign(phi(w_ind)) .* min(...
    -1 * [adj{1}(w_ind).*abs(s{1}(w_ind)), adj{2}(w_ind).*abs(s{2}(w_ind)), ...
    adj{3}(w_ind).*abs(s{3}(w_ind)), adj{4}(w_ind).*abs(s{4}(w_ind))], [], 2);
w = sparse(dims(1), dims(2));
capped_dp = min([abs(dp(w_ind)), (0.5-1e-1)*ones(size(w_ind))], [], 2);
w(w_ind) = -phi_i ./ (2*capped_dp - 1) - phi(w_ind);
% w(w_ind) = -phi_i - phi(w_ind);


    %
    % Solve the following problem:
    %     minimize ||dphi||^2
    %     subject to C dphi == dp
    % 

C = dp_dg * dg_dphi * S_phi';
d = dp(:);
b = zeros(size(C,2), 1) + S_phi * w(:);
% b = zeros(size(C,2), 1);
x = lso_priv_solver(b, C, d);
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
