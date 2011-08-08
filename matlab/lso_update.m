function [phi] = lso_update(phi, dp, p_eval, max_isles, step_sizes, varargin)
% PHI = LSO_UPDATE(PHI, DP, P_EVAL, MAX_ISLES, STEP_SIZES, [SEL])
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
%     MAX_ISLES: Non-negative integer.
%         Maximum number of islands to nucleate at each iteration.
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

dphi = lso_updatedir(phi, dp, sel);
[phi, dphi_isles] = lso_islands(phi, dp, max_isles, sel);
my_add = @(s) phi + s * (dphi + dphi_isles);


    %
    % Try various step sizes.
    %

step_sizes = sort(step_sizes(:), 'descend');
for k = 1 : length(step_sizes)
    if (p_eval(lso_fracfill(my_add(step_sizes(k)))) < eval0)
        % Form the updated level-set function.
        phi = my_add(step_sizes(k));
        return % If we have decreased the error, then stick with this step size.
    end
end

% If we get here, then we have not been able to decrease the error.
warning('Unable to decrease error function.');
