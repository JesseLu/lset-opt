function lso_test_updatedir(time_limit)
% LSO_TEST_UPDATEDIR(TIME_LIMIT)
% 
% Description
%     Run tests on LSO_UPDATEDIR for TIME_LIMIT seconds.

% Record results of tests.
dims = [];
err = [];

% Start tests.
start_time = tic;
while (toc(start_time) < time_limit)
    dims(end+1,:) = randi(10, [1 2]); % Pick dimensions of grid.

    % Choose variables.
    phi = lso_regularize(rand(dims(end,:)) - rand(1));
    p = lso_fracfill(phi);
    dp = randn(dims(end,:));
    dp = dp ./ norm(dp);

    % Get the update direction and compute error.
    dphi = lso_updatedir(phi, dp);
    err(end+1) = max(max(abs(((p > -1) & (p < 1)) .* ...
        (lso_fracfill(phi + 1e-5 * dphi) - (p + 1e-5 * dp)))));
end

% Print out results.
for k = 1 : length(err)
    fprintf('(%d, %d):\t%e\t', dims(k,:), err(k));
    if (err(k) < 1e-6)
        fprintf('passed.');
    else
        fprintf('FAILED!');
    end
    fprintf('\n');
end
