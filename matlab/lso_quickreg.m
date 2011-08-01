function [phi] = lso_quickreg(phi_hat)
% phi = lso_quickreg(phi)
% 
% Description
%     Quick regularization.

    % 
    % Step 1: Set all phi_hat == 0 to phi_hat = eps.
    %

phi_hat = phi_hat + eps * (phi_hat == 0);


    %
    % Step 2: Fix cells which are not adjacent to boundaries.
    %

[adj, on_border] = lso_priv_adjacents(phi_hat); % Find border cells.
phi = (~on_border) .* sign(phi_hat) + (on_border) .* phi_hat; % Fix values of non-border cells.

