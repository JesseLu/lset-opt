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
    dims(end+1,:) = randi(30, [1 2]); % Pick dimensions of grid.

    % Choose variables.
    phi = lso_regularize(rand(dims(end,:)) - rand(1));
end
