function [p] = lso_fracfill(phi)
% P = LSO_FRACFILL(PHI)
% 
% Description
%     Calculate the fractional-filling of the grid based on the level-set
%     function.
% 
% Input
%     PHI: 2d array (level-set function).
% 
% Output
%     P: 2d array (fractional filling).
%         The values of P range from -1 to 1.

    %
    % Calculate gammas.
    %

gamma = lso_priv_gamma(phi);


    % 
    % Calculate fractional-filling.
    %

p = 0.5 * sign(phi) .* ((gamma{1} + gamma{2} - 1) + (gamma{3} + gamma{4} - 1));

    

