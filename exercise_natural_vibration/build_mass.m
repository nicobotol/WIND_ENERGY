function [M] = build_mass(r_item)
% This funcition builds the mass matrix of the turbine


parameters;

% load data for the static deflaction
static_blade_data = readmatrix(structural_filename);
[r_vector, ~, ~, ~] = load_blade_data(blade_filename);

M = zeros(2*r_item);

% for i=2:r_item-1
%   l1 = 0.5*(r_vector(i) - r_vector(i - 1)); % length of blades at the left of the node
%   l2 = 0.5*(r_vector(i + 1) - r_vector(i)); % length of the blade at the right of the node
%   M(2*i - 3, 2*i - 3) = static_blade_data(i, 3)*(l1 + l2);
%   M(2*i - 2, 2*i - 2) = static_blade_data(i, 3)*(l1 + l2);
% end
% % mass of the last element
% M(2*r_item - 3, 2*r_item - 3) = 0.5*(r_vector(end) - r_vector(end - 1))*static_blade_data(end, 3);
% M(2*r_item - 2, 2*r_item - 2) = 0.5*(r_vector(end) - r_vector(end - 1))*static_blade_data(end, 3);

for i=1:r_item
  M(2*i - 1, 2*i - 1) = static_blade_data(i, 3);
  M(2*i, 2*i) = static_blade_data(i, 3);
end
% mass of the last element
% M(2*r_item - 3, 2*r_item - 3) = static_blade_data(end, 3);
% M(2*r_item - 2, 2*r_item - 2) = static_blade_data(end, 3);


end