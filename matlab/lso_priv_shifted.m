function [s] = lso_priv_shifted(a)
% S = LSO_PRIV_SHIFTED(A)
% 
% Description
%     Shift a 2d array to one-cell in each direction.
% 
% Inputs
%     A: 2d array.
% 
% Outputs
%     S: 4-element cell of 2d arrays.
%         S{1,2,3,4} are shifted arrays in the right, left, below, and above
%         directions. Zero-padding is used so that array-size is same as A.


dims = size(a);
s{1} = cat(1, a(2:end,:), zeros(1,dims(2))); % To the right.
s{2} = cat(1, zeros(1,dims(2)), a(1:end-1,:)); % To the left.
s{3} = cat(2, a(:,2:end), zeros(dims(1),1)); % Downwards.
s{4} = cat(2, zeros(dims(1),1), a(:,1:end-1)); % Upwards.
