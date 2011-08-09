function [phi] = lso_cheese(dims)
% PHI = LSO_CHEESE(DIMS)
% 
% Description
%     Create a checkerboard level-set pattern.
% 
% Inputs
%     DIMS: 2-element array of positive integers.
%         Size of grid.
% 
% Outputs
%     PHI: 2d array (level-set function).
        

% Create alternating pattern.
phi = repmat([1 -1; -1 1], ceil(dims(1)/2), ceil(dims(2)/2)); 
phi = phi(1:dims(1), 1:dims(2)); % Trim.
