function [phi] = lso_update(phi, dp, p_eval, step_sizes, varargin)
% PHI = LSO_UPDATE(PHI, DP, P_EVAL, STEP_SIZES, [SEL])
% 
% Description
%     Update the level-set function according to the derivative in the 
%     fractional-filling.
% 
% Inputs
%     PHI: 2d array (level-set function).
% 
%     DP: 2d array.
%         The derivative in the fractional-filling values.
% 
%     P_EVAL: Function handle.
%         P_EVAL must accept a 2D fractional-filling array (p) and return a 
%         non-negative scalar. LSO_UPDATE guarantees that the output PHI
%         produces less error than the input value of PHI.
% 
%     STEP_SIZES: Array of non-negative numbers.
%         The various step sizes with which to attempt to decrease P_EVAL. 
%         These do not need to be ordered.
%
%     SEL: 2d array (optional).
%         Marks the active (changeable) cells with a 1 and inactive cells with
%         0. Default value: all cells active.
% 
% Outputs
%     PHI: 2d array (level-set function).
%         Updated level-set function.


% Initial error, we need to decrease this to take a valid step.
eval0 = p_eval(lso_fracfill(phi));

% Determine the active (non-fixed) cells.
if isempty(varargin)
    sel = ones(size(phi));
else
    sel = varargin{1};
end
 

    %
    % Calculate dphi based on dp.
    %

dphi = zeros(size(phi));


    % 
    % Find cells which are next to a border and active.
    %

[adj, on_border] = lso_priv_adjacents(phi); % Find cells on border.
ind = find(on_border & sel); % Indices of active, on-border cells.

% If no active on-border cells, then we don't need to change phi.
if isempty(ind)
    return
end
 
    %
    % Calculate dphi using the gradient.
    %

grad = phi ./ lso_priv_gamma(phi);
dphi(ind) = dp(ind) .* abs(grad(ind));


    %
    % Try various step sizes.
    %

step_sizes = sort(step_sizes(:), 'descend'); % Try largest step first.

for k = 1 : length(step_sizes)
    if (p_eval(lso_fracfill(phi + step_sizes(k) * dphi)) < eval0)
        phi = phi + step_sizes(k) * dphi; % Form the updated level-set function.
        return 
    end
end

% If we get here, then we have not been able to decrease the error.
warning('Unable to decrease error function.');
