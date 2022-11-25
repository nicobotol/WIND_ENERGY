% function [py, pz, Theta_p] = internal_beam_loads_eigen(V0, n)
function [py, pz] = internal_beam_loads_eigen(n)
% This function computes the internal loads on the blades, modeling it as a
% beam

% load data from parameters
parameters

% load the blade data from "bladedat.txt"
[r_vector, ~, ~, ~] = load_blade_data(blade_filename);

r_item = size(r_vector, 2); % number of cross sections along the blade

% r_item_no_tip = r_item - 1;

% initialize a vector with 0 and the size of the total number of points
% except the 1st
py = zeros(r_item, 1);
pz = zeros(r_item, 1);

% write 1 in the dof loaded
if ~mod(n, 2) == 1 % if the gdl is even, load on z
  pz(n/2) = 1;
else % if the gdl is odd, load on y
  py(floor(n/2) + 1) = 1;
end

end