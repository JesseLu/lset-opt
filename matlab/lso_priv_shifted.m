function [s] = lso_priv_shifted(a)

dims = size(a);
s{1} = cat(1, a(2:end,:), zeros(1,dims(2))); % To the right.
s{2} = cat(1, zeros(1,dims(2)), a(1:end-1,:)); % To the left.
s{3} = cat(2, a(:,2:end), zeros(dims(1),1)); % Downwards.
s{4} = cat(2, zeros(dims(1),1), a(:,1:end-1)); % Upwards.
