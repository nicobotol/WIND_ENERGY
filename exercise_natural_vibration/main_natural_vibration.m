clear 
close all
clc

%% Load data
% load the physical parameters from the file "parameters.m"
parameters

% load pitchig angle
load('Pitching.mat'); % first row: windspeed between 4-25 (m/s)
                      % second row: angle for stall
                      % third row: angle for feathering

% load the airfoil data from the different files
[aoa_mat, cl_mat, cd_mat] = load_airfoil_data(filenames);

% load the blade data from "bladedat.txt"
[r_vector, c_vector, beta_vector, thick_vector] = load_blade_data(blade_filename);

r_item = size(r_vector, 2); % number of cross sections along the blade
r_item_no_tip = r_item - 1; % numer of cross section to take into account

% build the mass matrix
[M] = build_mass(r_item);

% build the flexibility matrix
F = zeros(2*r_item);
for n=1:2*r_item
 [u] = deflection_eigen(r_item, n);  
 F(:, n) = u;
 
end

% solve the eigenvalue problem
[V, D] = eigs(F*M);

% find the eigenvalues
omega1 = sqrt(1/D(1,1));
omega2 = sqrt(1/D(2,2));
omega3 = sqrt(1/D(3,3));
omega4 = sqrt(1/D(4,4));
omega5 = sqrt(1/D(5,5));
omega6 = sqrt(1/D(6,6));
f1 = omega1/(2*pi);
f2 = omega2/(2*pi);
f3 = omega3/(2*pi);

fprintf("Mode 1: %f [rad/s] - %f [Hz] \n", omega1, f1); 
fprintf("Mode 2: %f [rad/s] - %f [Hz] \n", omega2, f2); 
fprintf("Mode 3: %f [rad/s] - %f [Hz] \n \n", omega3, f3); 
% find the eigenvectors
MS1_y = V(1:2:end, 1)/max(abs(V(:, 1)));
MS1_z = V(2:2:end, 1)/max(abs(V(:, 1)));

MS2_y = V(1:2:end, 2)/max(abs(V(:, 2)));
MS2_z = V(2:2:end, 2)/max(abs(V(:, 2)));

MS3_y = V(1:2:end, 3)/max(abs(V(:, 3)));
MS3_z = V(2:2:end, 3)/max(abs(V(:, 3)));

%% 
fig_eigenmodes = figure('Position', get(0, 'Screensize'));
subplot(3,1,1)
plot(r_vector, MS1_y, 'LineWidth',line_width )
hold on
plot(r_vector, MS1_z,'LineWidth',line_width )
hold off
ylabel('u/u_{MAX} [-]')
grid on
legend('u_y', 'u_z', 'Location', 'southwest')
title(strcat('1^{st} flap - \omega = ', num2str(omega1), ' [rad/s]', ...
  ' f = ', num2str(f1), ' [Hz]') )
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)

subplot(3,1,2)
plot(r_vector, MS2_y, 'LineWidth',line_width )
hold on
plot(r_vector, MS2_z, 'LineWidth',line_width )
hold off
ylabel('u/u_{MAX} [-]')
grid on
title(strcat('1^{st} edge - \omega = ', num2str(omega2), ' [rad/s]', ...
  ' f = ', num2str(f2), ' [Hz]') )
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)

subplot(3,1,3)
plot(r_vector, MS3_y,'LineWidth',line_width )
hold on
plot(r_vector, MS3_z, 'LineWidth',line_width )
hold off
xlabel('r [m]')
ylabel('u/u_{MAX} [-]')
grid on
title(strcat('2^{nd} flap - \omega = ', num2str(omega3), ' [rad/s]', ...
  ' f = ', num2str(f3), ' [Hz]') )
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_eigenmodes, ['C:\Users\Niccolò\Documents\UNIVERSITA\' ...
  '5° ANNO\WIND_ENERGY\exercise_natural_vibration\' ...
  'fig_eigenmodes.png'],'png');

w1 = 3.516*(1/1)^0.5/(18-1)^2;
w2 = 22.03*(1/1)^0.5/(18-1)^2;
w3 = 61.7*(1/1)^0.5/(18-1)^2;

%% Trust force
if V0 < V0_rated - 1e-3 % velocities below the rated one
  lambda = lambda_opt;
  omega = V0 * lambda / R;
  Theta_p = Theta_p_opt; % pitch    
else
  lambda = omega_max * R / V0;
  omega = omega_max;
  Theta_p = interp1(Theta_p_limit(1,:), Theta_p_limit(3,:), V0); % pitch
end

[cp_partial, cT_partial, pt, pn] = cP_cT_partial(r_item_no_tip, r_vector, ...
beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, R, lambda, ...
Theta_p, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0, ...
omega, i_max);
     
cP = lambda*B/(R*A)*trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial);
cT = B/A*trapezoidal_integral(r_vector(1:r_item_no_tip), cT_partial);
      
P = cP*0.5*A*rho*V0^3;
T = cT*0.5*A*rho*V0^2/1000;