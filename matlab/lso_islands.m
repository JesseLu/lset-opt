function [ind] = lso_islands(phi, dp)

dims = size(phi);


    %
    % Find which cells are viable for island nucleation.
    %

% Shifted fractional-filling values.
p = lso_fracfill(phi);
s{1} = cat(1, p(2:end,:), p(end,:)); % To the right.
s{2} = cat(1, p(1,:), p(1:end-1,:)); % To the left.
s{3} = cat(2, p(:,2:end), p(:,end)); % Downwards.
s{4} = cat(2, p(:,1), p(:,1:end-1)); % Upwards..

% Need the current and four adjacent cells to be all either 1 or -1.
viable = ...
    ((dp < 0) & (p==1) & (s{1}==1) & (s{2}==1) & (s{3}==1) & (s{4}==1)) | ...
    ((dp > 0) & (p==-1) & (s{1}==-1) & (s{2}==-1) & (s{3}==-1) & (s{4}==-1));


    %
    % Go through and mark as many cells as possible for nucleation,
    % sort by strength of dp at that cell.
    %

ind = [];
while(sum(viable(:)) > 0)
% for cnt = 1 : 2
    % Find viable cell with strongest dp.
    [temp, ind(end+1)] = max(viable(:) .* abs(dp(:)));

    % Eliminate neighboring viable cells.
    [i, j] = ind2sub(dims, ind(end));
    viable(array_access(dims, i + [0 -1 1 0 0], j + [0 0 0 -1 1])) = 0;
end


function [ind] = array_access(dims, i, j)
cap  = @(k, len) k .* ((k >= 1) & (k <= len)) + 1 * (k < 1) + len * (k > len);
ind = sub2ind(dims, cap(i, dims(1)), cap(j, dims(2)));
