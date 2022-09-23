function [pn_ashes, pt_ashes, r_ashes] = import_ashes()
% Import data from the ashes code

pn_ashes = zeros(4, 40);
pt_ashes = zeros(4, 40);
r_ashes = zeros(1, 40);

pn_load = load("ashes_data\p_n.mat");
pn_ashes = pn_load.p_n_matrix(:, :);

pt_load = load("ashes_data\p_t.mat");
pt_ashes = pt_load.p_t_matrix(:, :);

bladepos = load("ashes_data\bladepos.mat");
r_ashes = bladepos.blade_pos(1,:) + 2.8; % load the radial position

end