%% This code computes the BEM code for one entire blade, and not only for a section
% Inputs parameters for the code to work, to be setted in parameters.m are
% the wind speed V0, the rotational speed omega, and the pitch angle 

clear 
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMETERS IN parameters.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load the physical parameters from the file "parameters.m"
parameters

% load the airfoil data from the different files
[aoa_mat, cl_mat, cd_mat] = load_airfoil_data(filenames);

% load the blade data from "bladedat.txt"
[r_vector, c_vector, beta_vector, thick_vector] = ...
load_blade_data(blade_filename);

r_item = size(r_vector, 2); % number of cross sections along the blade
r_item_no_tip = r_item - 1; % numer of cross section to take into account

pitch_vector = linspace(pitch_range(1), pitch_range(2), pitch_item); % 
% vector of pitch equally distributed in the range 

lambda = omega * R / V0;

[cp_partial, ct_partial, pt, pn] = cP_cT_partial(r_item_no_tip, ...
  r_vector, beta_vector, thick_vector, c_vector, B, a_guess, ...
  a_prime_guess, R, lambda, Theta_p, aoa_mat, cl_mat, cd_mat, ...
  thick_prof, fake_zero, rho,V0,omega, i_max);

% compute the power
P = 0.5*omega*B*rho*V0^2*trapezoidal_integral( ...
  r_vector(1:r_item_no_tip), cp_partial);
T = 0.5*rho*V0^2*B*trapezoidal_integral( ...
  r_vector(1:r_item_no_tip), ct_partial);

fprintf("The power is: %f [MW]\n", P*1e-6);
fprintf("The thrust is %f [KN]\n", T*1e-3);


pt(end + 1) = 0;
pn(end + 1) = 0;

fig_loads = figure('Position', get(0, 'Screensize'));
plot(r_vector, pt, 'LineWidth',line_size)
hold on
plot(r_vector, pn, 'LineWidth',line_size)
hold off
grid on
legend('pt', 'pn', 'Location','northwest')
xlabel('Blade span [m]')
ylabel('Loads [N/m]')
title('Loads along the blade')
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_loads, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_BEM_entire_blade\fig_loads.png'],'png');
