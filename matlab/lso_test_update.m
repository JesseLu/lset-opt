function lso_test_update(time_limit)
% Tests for the lso_update function.

dims = [];
err = [];
start_time = tic;
while (toc(start_time) < time_limit)
    dims(end+1,:) = randi(100, [1 2]);
    phi = lso_regularize(rand(dims(end,:)) - rand(1));
    p = lso_fracfill(phi);
    dp = randn(dims(end,:));
    dp = dp ./ norm(dp);

    dphi = lso_update(phi, dp);
    err(end+1) = prod(dims(end,:)).^-1 * norm(((p > -1) & (p < 1)) .* ...
        (lso_fracfill(phi + 1e-5 * dphi) - (p + 1e-5 * dp)));
end

for k = 1 : length(err)
    fprintf('(%d, %d):\t%e\n', dims(k,:), err(k));
end
