function [err] = lso_test_updatedir(time_limit)
% LSO_TEST_UPDATEDIR(TIME_LIMIT)
% 
% Description
%     Run tests on LSO_UPDATEDIR for TIME_LIMIT seconds.

fail_lim = 0.01; % If error exceeds this, issue fail.  

% Start tests.
start_time = tic;

err = [];
while (toc(start_time) < time_limit)
    dims = randi(100, [1 2]); % Pick dimensions of grid.

    % Choose variables.
    phi = lso_regularize(rand(dims) - rand(1));
    p = lso_fracfill(phi);
    dp = randn(dims);
    dp = dp ./ norm(dp(:));

%     p1 = p + dp
%     dp = 

    % Get the update direction and compute error.
    dphi = lso_updatedir(phi, dp);
    s = 1e-5;
    e = abs(((p > -1) & (p < 1)) .* ...
        (lso_fracfill(phi + s * dphi) - (p + s * dp)) ./ s);
    err(end+1) = max(e(:)); 

    % Print results.
    fprintf('(%d, %d):\t%e\t', dims, err(end));

    if (err < fail_lim)
        fprintf('passed.\n');
    else
        [temp, ind] = max(e(:));
        [i, j] = ind2sub(size(e), ind);
        k = 2;
        subplot 121; imagesc(e(i-k:i+k,j-k:j+k)'); axis equal tight;
        f = phi + 1e-5 * dphi;
        subplot 122; lso_plot(f(i-k:i+k, j-k:j+k));
        phi(i,j)
        dphi(i,j)
        p(i,j)
        dp(i,j)

        fprintf('FAILED! error (%e)\n', max(err(:)));
        return
    end
end




