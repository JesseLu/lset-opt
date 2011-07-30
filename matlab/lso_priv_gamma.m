function [gamma] = lso_priv_gamma(phi)

dims = size(phi); % Size of the array.


    %
    % Calculate gammas.
    %

% Shifted arrays of phi.
phi_shift = {phi(1:end-1,:), phi(2:end,:), phi(:,1:end-1), phi(:,2:end)};

% Calculate the values of the gammas.
gamma{1} = cat(1, ... % To the right.
    my_gamma(phi_shift{1}, phi_shift{2}) , ...
    ones(1, dims(2)));

gamma{2} = cat(1, ... % To the left.
    ones(1, dims(2)), ...
    my_gamma(phi_shift{2}, phi_shift{1}));

gamma{3} = cat(2, ... % Downward.
    my_gamma(phi_shift{3}, phi_shift{4}), ...
    ones(dims(1), 1));

gamma{4} = cat(2, ... % Upward.
    ones(dims(1), 1), ...
    my_gamma(phi_shift{4}, phi_shift{3}));


    %
    % Helps to calculate gamma.
    %

function [gamma] = my_gamma(phi_0, phi_i)
gamma = phi_0 ./ (phi_0 - phi_i); % Calculate distance.
s = sign(phi_0) ~= sign(phi_i); % Eligible or not?
gamma(find(~s)) = 1.0; % If ineligible, replace with distance of 1.

