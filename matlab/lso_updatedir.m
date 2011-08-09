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
dphi = zeros(dims);


    % 
    % Find cells which are next to a border and active.
    %

[adj, on_border] = lso_priv_adjacents(phi); % Find cells on border.
ind = find(on_border & sel); % Indices of active, on-border cells.

% If no active on-border cells, return dphi = 0.
if isempty(ind)
    return
end
 
    %
    % Calculate dphi using the gradient.
    %

grad = phi ./ lso_priv_gamma(phi);
dphi(ind) = dp(ind) .* abs(grad(ind));


