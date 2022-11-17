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



%% Internal loads
[py, pz, Theta_p] = internal_beam_loads(V0);

fig_load = figure('Position', get(0, 'Screensize'));
plot(r_vector, pz, 'LineWidth', line_width)
hold on
plot(r_vector, py, 'LineWidth', line_width)
hold off
legend('pz', 'py', 'Location','northwest')
xlabel('Spanwise length [m]')
ylabel('Load [N]')
grid on
set(gca, 'FontSize', font_size);
saveas(fig_load, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\exercise_static_deflaction\figures\fig_load.png','png');

%% Deflaction and angle

[Theta_y, Theta_z, u_y, u_z, Ky, Kz] = deflection(V0, r_item);

fig_angle = figure('Position', get(0, 'Screensize'));
plot(r_vector, rad2deg(Theta_y), 'LineWidth', line_width)
hold on
plot(r_vector, rad2deg(Theta_z), 'LineWidth', line_width)
hold off
legend('\theta_y', '\theta_z', 'Location','southwest')
xlabel('Spanwise position [m]')
ylabel('Deflaction angle [°]')
grid on
set(gca, 'FontSize', font_size);
saveas(fig_angle, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\exercise_static_deflaction\figures\fig_angle.png','png');

fig_deflection = figure('Position', get(0, 'Screensize'));
plot(r_vector, u_y, 'LineWidth', line_width)
hold on
plot(r_vector, u_z, 'LineWidth', line_width)
hold off
legend('u_y', 'u_z', 'Location','northwest')
xlabel('Spanwise position [m]')
ylabel('Deflection [m]')
grid on
set(gca, 'FontSize', font_size);
saveas(fig_deflection, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\exercise_static_deflaction\figures\fig_deflection.png','png');

fig_K = figure('Position', get(0, 'Screensize'));
plot(r_vector, rad2deg(Ky), 'LineWidth', line_width)
hold on
plot(r_vector, rad2deg(Kz), 'LineWidth', line_width)
hold off
legend('K_y', 'K_z', 'Location','southwest')
xlabel('Spanwise position [m]')
ylabel('K [°/m]')
grid on
set(gca, 'FontSize', font_size);
saveas(fig_K, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\exercise_static_deflaction\figures\fig_K.png','png');