function [phi] = lso_update(phi, dp)

    %
    % Change phi relative to dp.
    %

dphi = lso_updatedir(phi, dp);
phi = phi + 3e-1 * dphi;


    %
    % Nucleate and form islands.
    %

phi = lso_islands(phi, dp);

    %
    % Check errors? Regularize?
    %
