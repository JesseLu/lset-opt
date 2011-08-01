function [phi] = lso_quickreg(phi_hat)
% PHI = LSO_QUICKREG(PHI)
% 
% Description
%     Perform a quick regularization of the proposed level-set function. In 
%     contrast to LSO_REGULARIZE, no matrix-solve is required.
% 
%     A quick regularization involves only the first two steps of normal
%     regularization. That is, any PHI_HAT = 0 are set to PHI = EPS, where 
%     EPS is the lowest positive number available. Also, non-border values are
%     set to either +1 or -1.
% 
% Inputs
%     PHI_HAT: 2d array (level-set function).
% 
% Outputs
%     PHI: 2d array (level-set function).


    % 
    % Step 1: Set all phi_hat == 0 to phi_hat = eps.
    %

phi_hat = phi_hat + eps * (phi_hat == 0);


    %
    % Step 2: Fix cells which are not adjacent to boundaries.
    %

[adj, on_border] = lso_priv_adjacents(phi_hat); % Find border cells.
phi = (~on_border) .* sign(phi_hat) + (on_border) .* phi_hat; 

