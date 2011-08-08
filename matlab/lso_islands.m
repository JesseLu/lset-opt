function [phi, dphi] = lso_islands(phi, dp, max_isles, sel)
% [PHI, DPHI] = LSO_ISLANDS(PHI, DP, MAX_ISLES, SEL)
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
%     MAX_ISLES: non-negative integer.
%         The maximum number of new islands to nucleate.
%
%     SEL: 2d array.
%         Marks the active (changeable) cells with a 1 and inactive cells with
%         0. 
%             
% Outputs
%     PHI: 2d array (level-set function).
%         New level-set function which is set up for island nucleation. 
%         However, no islands are nucleated on PHI yet.
% 
%     DPHI: 2d array.
%         PHI + a * DPHI, where a > 0, results in a level-set function with 
%         nucleated islands.

    
    %
    % Check inputs.
    %

dims = size(phi);
dphi = zeros(dims); % Output variable.

if (max_isles <= 0) % If we don't want any islands, we're done!
    return
end


    %
    % Find which cells are viable for island nucleation.
    %

s = lso_priv_shifted(phi);

% Need the values of phi of the current and four adjacent cells to be all 
% either 1 or -1. Also, the value of dp must be of the opposite sign.
viable = sel &  ... 
    ((dp < 0) & (phi==1) & (s{1}==1) & (s{2}==1) & (s{3}==1) & (s{4}==1)) | ...
    ((dp > 0) & (phi==-1) & (s{1}==-1) & (s{2}==-1) & (s{3}==-1) & (s{4}==-1));

if isempty(find(viable)) % If no viable spots, exit without changing phi.
    return
end


    %
    % Go through and mark as many cells as possible for nucleation,
    % sort by strength of dp at that cell.
    %

% Calculate the strength of dp at the current cell.
% Uses a heuristic to approximate smearing.
t = lso_priv_shifted(dp); 
dp_stren = norm([4 1 1 1 1])^-1 * (4*dp + t{1} + t{2} + t{3} + t{4});

ind = []; % Indices of phi which are viable.
num_isles = 0;
while (sum(viable(:)) > 0) & (num_isles < max_isles)
    % Find viable cell with strongest dp_stren.
    [temp, ind(end+1)] = max(viable(:) .* abs(dp_stren(:)));

    % Make sure dp_stren is pointing in the same direction as dp
    if sign(dp(ind(end))) ~= sign(dp_stren(ind(end)))
        break % If it is not, then we're done finding islands to nucleate.
    end

    % Eliminate neighboring viable cells.
    viable(ind(end) + [0 -1 1 -dims(1) dims(1)]) = 0;
    
    num_isles = num_isles + 1; % Record that we nucleated one island.
end


    %
    % Form islands to exactly match capped dp values.
    %

% Cap relevant dp values at 1.
dp_cap = abs(dp(ind));
dp_cap = 1 * (dp_cap >= 1) + dp_cap .* (dp_cap < 1);

% How to change phi to form islands.
dphi(ind) = (dp_cap ./ (dp_cap - 2)) .* sign(phi(ind));

% New baseline phi, ready for island nucleation.
phi(ind) = eps .* sign(phi(ind)); 

