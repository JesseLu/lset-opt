function [dphi] = lso_update(phi, dp)

dims = size(phi);
N = prod(dims);


%     %
%     % Form selection matrix for p that are next to a border.
%     % These are the cells whose fractional-filling can actually be perturbed.
%     %
% 
gamma = lso_priv_gamma(phi); % Find boundary points.
% ind = find((gamma{1} ~= 0.5) | (gamma{2} ~= 0.5) | ... % Cells not on boundary.
%     (gamma{3} ~= 0.5) | (gamma{4} ~= 0.5));
% S_p = sparse(1:length(ind), ind, ones(length(ind), 1), length(ind), N);


    %
    % Form the dp / dgamma matrix.
    %

my_diag = @(x) spdiags(x(:), 0, numel(x), numel(x));
width = 2 * sign(phi) .* (gamma{1} + gamma{2});
height = 2 * sign(phi) .* (gamma{3} + gamma{4});
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
    %     subject to A dphi == dp
    % 

C = dp_dg * dg_dphi * S_phi';
d = dp(:);
b = zeros(size(C,2), 1);
x = b - C' * ((C * C') \ (C * b - d));
dphi = reshape(S_phi' * x, size(phi));





function [A] = my_dg_dphi(phi, adj, s)
N = numel(phi);
ind = find(adj);
if (s > 0)
    phi_ind = ind + s * (abs(phi(ind)) >= abs(phi(ind+s)));
else
    phi_ind = ind + s * (abs(phi(ind)) > abs(phi(ind+s)));
end
A = sparse(ind, phi_ind, -phi(phi_ind)./(phi(ind)-phi(ind+s)).^2, N, N);
