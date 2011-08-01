function [phi] = lso_update(phi, sel_phi, dp, p_eval, num_isles)

eval0 = p_eval(lso_fracfill(phi));


    %
    % Change phi relative to dp.
    %

dphi = lso_updatedir(phi, dp);
[phi, dphi_isles] = lso_islands(phi, dp, num_isles);
dphi = dphi + dphi_isles;
dphi = sel_phi .* dphi;
s = 2^5;
while p_eval(lso_fracfill(phi + s * dphi)) > eval0
    s = s / 2;
end

phi = phi + s * dphi;

% phi = lso_islands(phi, sel_phi .* dp);


    %
    % Check errors? Regularize?
    %
