function lso_test_updatedir(time_limit)
% LSO_TEST_UPDATEDIR(TIME_LIMIT)
% 
% Description
%     Run tests on LSO_UPDATEDIR for TIME_LIMIT seconds.

fail_lim = 1e-6; % If error exceeds this, issue fail.  

% Start tests.
start_time = tic;

while (toc(start_time) < time_limit)
    dims = randi(100, [1 2]); % Pick dimensions of grid.

    % Choose variables.
    phi = lso_regularize(rand(dims) - rand(1));
    p = lso_fracfill(phi);
    dp = randn(dims);
    dp = dp ./ norm(dp);

    % Get the update direction and compute error.
    dphi = lso_updatedir(phi, dp);
    err = abs(((p > -1) & (p < 1)) .* ...
        (lso_fracfill(phi + 1e-5 * dphi) - (p + 1e-5 * dp)));

    % Print results.
    fprintf('(%d, %d):\t', dims);
    if (max(err(:)) < fail_lim)
        fprintf('passed.\n');
    else
        subplot 121; imagesc(err'); axis equal tight;
        subplot 122; lso_plot(phi + 1e-5 * dphi);

        fprintf('FAILED! error (%e)\n', max(err(:)));
        return
    end
end

