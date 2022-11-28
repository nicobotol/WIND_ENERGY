function [py, pz, Theta_p] = internal_beam_loads(V0)
% This function computes the internal loads on the blades, modeling it as a
% beam

% load data from parameters
parameters

% load pitchig angle
load('Pitching.mat'); % first row: windspeed between 4-25 (m/s)
                      % second row: angle for stall
                      % third row: angle for feathering


% load the blade data from "bladedat.txt"
[r_vector, c_vector, beta_vector, thick_vector] = load_blade_data(blade_filename);

% load the airfoil data from the different files
[aoa_mat, cl_mat, cd_mat] = load_airfoil_data(filenames);

r_item = size(r_vector, 2); % number of cross sections along the blade

r_item_no_tip = r_item - 1;

pz = zeros(r_item, 1);
py = zeros(r_item, 1);
  
if V0 < V0_rated - 1e-3 % velocities below the rated one
  lambda = lambda_opt;
  omega = V0 * lambda / R;
  Theta_p = Theta_p_opt;      

else
  lambda = omega_max * R / V0;
  omega = omega_max;
  
  % pitch
  Theta_p = interp1(Theta_p_limit(1,:), Theta_p_limit(3,:), V0);

end

[~, ~, py, pz] = cP_cT_partial(r_item_no_tip, r_vector, ...
  beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, ...
  R, lambda, Theta_p, aoa_mat, cl_mat, cd_mat, thick_prof, ...
  fake_zero, rho, V0, omega, i_max);
pz(end + 1) = 0; % add value correspnding to the tip
py(end + 1) = 0;  % add value correspnding to the tip

end