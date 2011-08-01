function lso_test_regularize(time_limit)
% LSO_TEST_REGULARIZE(TIME_LIMIT)
% 
% Description
%     Run tests on LSO_REGULARIZE for TIME_LIMIT seconds.
% 
% Example
%     lso_test_regularize(10); % Test lso_regularize for 10 seconds.

fail_lim = 1e-6; % If error is more than this, issue a fail.

% Start tests.
start_time = tic;
while (toc(start_time) < time_limit)
    dims = randi(100, [1 2]); % Pick dimensions of grid.
    fprintf('(%d, %d): ', dims);

    % Choose variables.
    [phi, err_bnd, err_p] = lso_regularize(rand(dims) - rand(1));
    if (err_bnd < fail_lim) & (err_p < fail_lim)
        fprintf('passed.\n');
    else
        fprintf('FAILED!\n');
        return
    end
end
