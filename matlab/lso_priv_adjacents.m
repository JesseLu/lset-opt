function [adj, on_border] = lso_priv_adjacents(phi)
% [ADJ, ON_BORDER] = LSO_PRIV_ADJACENTS(PHI)
% 
% Description
%     Finds whether a cell is adjacent to a border on either of it's four
%     borders.
% 
% Inputs
%     PHI: 2d array (level-set function).
% 
% Outputs
%     ADJ: 4-element cell of 2d arrays.
%     ADJ{I} for I = 1,2,3,4 denote whether there is a border to the right,
%     left, below, and above the current cell respectively. A value of 1 means
%     that a border is present, while a value of 0 means that one is not.
% 
%     ON_BORDER: 2d array.
%     1 if cell is adjacent to a border (in some direction), 0 otherwise.


dims = size(phi); % Dimensions of the grid.

% Boundary present to the right and left of the cell
horiz = sign(phi(1:end-1,:)) ~= sign(phi(2:end,:));
adj{1} = cat(1, horiz, zeros(1, dims(2))); % To the right.
adj{2} = cat(1, zeros(1, dims(2)), horiz); % To the left.

% Boundary present below and above the cell
vert = sign(phi(:,1:end-1)) ~= sign(phi(:,2:end));
adj{3} = cat(2, vert, zeros(dims(1), 1)); % Downwards.
adj{4} = cat(2, zeros(dims(1), 1), vert); % Updwards.

% Find which cells are next to a boundary point.
on_border = adj{1} | adj{2} | adj{3} | adj{4};
