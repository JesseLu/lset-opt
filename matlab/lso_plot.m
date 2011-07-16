function lso_plot(phi)
% LSO_PLOT(PHI)
% 
% Description
%     Simply plots the filling fraction, as well as the boundary points of the
%     level-set function.
% 
% Input
%     PHI: 2d array (level-set function).
% 
% Output
%     None.
% 
% Example
%     phi_hat = rand(40, 30) - 0.2;
%     lso_plot(phi_hat)   


% Plot the filling fraction.
imagesc(lso_fracfill(phi)', [-1 1]); 
a = axis; 

% Plot the boundary points.
hold on
[x, y] = lso_boundaries(phi);
plot(x, y, 'r.');
hold off

% Make things look pretty.
axis(a);
axis equal tight;
