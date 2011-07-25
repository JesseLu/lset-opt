function [phi] = lso_update(phi, dp)

    %
    % Change phi relative to dp.
    %

dphi = lso_updatedir(phi, dp);
phi = phi + dphi;


    %
    % Nucleate and form islands.
    %

isles = lso_islands(phi, dp);

    %
    % Check errors? Regularize?
    %
