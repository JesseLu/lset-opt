function lso_test_update(num_isles)
% 
% Description
%     Run tests on LSO_ISLANDS for TIME_LIMIT seconds.


dims = randi(10, [1 2])+10; % Pick dimensions of grid.
fprintf('(%d, %d):\n', dims);


phi_target = randn(dims);
d = 4;
phi_target = interp2(phi_target, [1:d^-1:dims(2)]', 1:d^-1:dims(1));
dims = size(phi_target)

p_target = lso_fracfill(phi_target);
sel_phi = zeros(dims);
sel_phi(2:end-1, 2:end-1) = 1;

phi = (1-sel_phi) .* (phi_target) - sel_phi;
dp = @(p) (p_target - p);
err = @(p) norm(p_target(:) - p(:));

for k = 1 :1e4 
    
    subplot 121; lso_plot(phi); 
    subplot 122; lso_plot(phi_target); 
    % subplot 122; imagesc((p_target - lso_fracfill(phi))'); 
    axis equal tight;

    pause(0.4);
    phi = lso_update(phi, ones(dims), dp(lso_fracfill(phi)), err, num_isles);
    e = err(lso_fracfill(phi))
end


return

% Start tests.
start_time = tic;
while (toc(start_time) < time_limit)
    dims = randi(100, [1 2]); % Pick dimensions of grid.
    fprintf('(%d, %d):  \t', dims);

    % Choose variables.
    pol = 2 * (randi(2,1) - 1.5);
    phi0 = lso_regularize(randn(dims) + pol);
    dp = 2*rand(dims) - 1;

    phi1 = lso_islands(phi0, dp); % Make islands.

    % Calculate error.
    p = {lso_fracfill(phi0), lso_fracfill(phi1)};
    e = (p{2} ~= p{1}) .* ((p{2} - p{1}) - dp);
    err = max(abs(e(:)));
    fprintf('%e\t', err);

    % Tell us whether we passed or not.
    if (err < fail_lim)
        fprintf('passed.\n');
    else
        subplot 121; lso_plot(phi0);
        subplot 122; lso_plot(phi1);
        fprintf('FAILED!\n');
        return
    end
end
