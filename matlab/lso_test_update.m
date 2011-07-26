function lso_test_update()
% 
% Description
%     Run tests on LSO_ISLANDS for TIME_LIMIT seconds.

fail_lim = 1e-15; % If error is more than this, issue a fail.

dims = randi(3, [1 2])+5; % Pick dimensions of grid.
fprintf('(%d, %d):\n', dims);


phi_target = randn(dims);
p_target = lso_fracfill(phi_target);
phi = lso_regularize(randn(dims)-1e3);
dp = @(p) (p_target - p)./1e0;
err = @(p) norm(p_target(:) - p(:));

for k = 1 :1e4 
    phi = lso_regularize(lso_update(phi, dp(lso_fracfill(phi))));
    e = err(lso_fracfill(phi))
    
    subplot 121; lso_plot(phi); 
    subplot 122; lso_plot(phi_target); 
    % subplot 122; imagesc((p_target - lso_fracfill(phi))'); 
    axis equal tight;

    pause(0.4);
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
