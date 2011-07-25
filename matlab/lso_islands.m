function [phi] = lso_islands(phi, dp)
% [PHI] = LSO_ISLANDS(PHI, DP)
% 
% Description
%     Form islands at viable cells in the grid. A viable cell is not adjacent
%     to any boundary points, neither are its adjacent cells.
% 
% Inputs
%     PHI: 2d array (level-set function).
% 
%     DP: 2d array.
%         The derivative in the fractional-filling (p).
%             
% Outputs
%     PHI: 2d array (level-set function).
%         New level-set function (with islands).


dims = size(phi);


    %
    % Find which cells are viable for island nucleation.
    %

s{1} = cat(1, phi(2:end,:), zeros(1,dims(2))); % To the right.
s{2} = cat(1, zeros(1,dims(2)), phi(1:end-1,:)); % To the left.
s{3} = cat(2, phi(:,2:end), zeros(dims(1),1)); % Downwards.
s{4} = cat(2, zeros(dims(1),1), phi(:,1:end-1)); % Upwards.

% Need the current and four adjacent cells to be all either 1 or -1, and the
% value of dp must be of the opposite sign.
viable = ...
    ((dp < 0) & (phi==1) & (s{1}==1) & (s{2}==1) & (s{3}==1) & (s{4}==1)) | ...
    ((dp > 0) & (phi==-1) & (s{1}==-1) & (s{2}==-1) & (s{3}==-1) & (s{4}==-1));

if isempty(find(viable)) % If no viable spots, exit without changing phi.
    return
end


    %
    % Go through and mark as many cells as possible for nucleation,
    % sort by strength of dp at that cell.
    %

ind = [];
while(sum(viable(:)) > 0)
    % Find viable cell with strongest dp.
    [temp, ind(end+1)] = max(viable(:) .* abs(dp(:)));

    % Eliminate neighboring viable cells.
    viable(ind(end) + [0 -1 1 -dims(1) dims(1)]) = 0;
end


    %
    % Form islands to exactly match capped dp values.
    %

% Cap relevant dp values at 2.
dp_cap = abs(dp(ind));
dp_cap = 2 * (dp_cap >= 2) + dp_cap .* (dp_cap < 2);

% Form islands.
phi(ind) = (dp_cap ./ (dp_cap - 4)) .* sign(phi(ind));

