function lso_test_islands(time_limit)
% LSO_TEST_ISLANDS(TIME_LIMIT)
% 
% Description
%     Run tests on LSO_ISLANDS for TIME_LIMIT seconds.
% 
% Example
%     lso_test_islands(1); % Run tests on lso_islands for 1 second.

fail_lim = 1e-15; % If error is more than this, issue a fail.

% Start tests.
start_time = tic;
while (toc(start_time) < time_limit)
    dims = randi(5, [1 2]) + 5; % Pick dimensions of grid.
    fprintf('(%d, %d):  \t', dims);

    % Choose variables.
    pol = 2 * (randi(2,1) - 1.5);
    phi0 = lso_regularize(randn(dims) + pol);
    dp = 2*rand(dims) - 1;

    [phi1, dphi] = lso_islands(phi0, dp); % Make islands.
    phi1 = phi1 + dphi;

    % Calculate error.
    % Error is only calculated for cells on which islands were nucleated.
    p = {lso_fracfill(phi0), lso_fracfill(phi1)};
    e = (phi0 ~= phi1) .* ((p{2} - p{1}) - dp);
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
