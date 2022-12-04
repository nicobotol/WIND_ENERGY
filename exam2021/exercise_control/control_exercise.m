clear 
close all
clc

%% CONSTANT PARAMETERS
time_max = 50; % maximum time to be studied (s)
delta_t = 0.005; % discretization time (s)
i_max = time_max / delta_t; % number of iterations

sincronous_velocity = 5; % (rad/s)
MG_coefficient = 2.55e5; % generator curve coefficient

R = 9.5; % rotor diameter (m)
rho = 1.225; % air density (kg/m^3)
V0 = 10; % constant wind speed (m/s)
V0_f = 10; % wind speed for free running turbine (m/s)
omega_initial = 1; % initial angular velocity (rad/s)
omega_initial_f = 1; % initial angular velocity free running turbine (rad/s)
I = 1.5e4; % rotor inertia (kgm^2)


line_width = 2;
font_size = 25; 
%% Q1 & Q2: PLOT OF THE LOADED DATA
% load the data form the table
cP_vs_lambda = load('Cp_Lambda.txt');
V0_vs_time = load('usim.dat');

figure()
plot(cP_vs_lambda(:,1), cP_vs_lambda(:,2), 'LineWidth', line_width)
xlabel('\lambda')
ylabel('cP')
title('cP as function of \lambda')
grid on

figure()
plot(V0_vs_time(:,1), V0_vs_time(:,2), 'LineWidth', line_width)
title('Turbolent wind')
xlabel('Time (s)')
ylabel('Wind velocity (m/s)')
grid on

%% Q4: COMPUTATION OF THE ANGULAR VELOCITY, ROTOR AND GENERATOR TORQUE 

% bild a vector with the constant velocity: it is needed to use
% control_function
V0_q1 = ones(1, i_max)*V0;
% compute the torque and the velocity
[Mr, Mg, omega, delta_t_vector] = control_function(i_max, ...
  omega_initial, R, V0_q1, rho, cP_vs_lambda, MG_coefficient, ...
  sincronous_velocity, delta_t, I);

% plot the time evolution of the roataional speed
fig_omega = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector, omega, 'LineWidth', line_width)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
title('Rotor angular velocity')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_omega, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_omega.png'],'png');

% plot the time evolution of the generator torque
fig_Mg = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector, Mg, 'LineWidth', line_width)
xlabel('Time [s]')
ylabel('Generator torque [Nm]')
title('Generator torque')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_Mg, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_Mg.png'],'png');

% plot the time evolution of the wind torque
fig_Mr = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector, Mr, 'LineWidth', line_width)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_Mr, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_Mr.png'],'png');

% plot generator and wind torque on the same plot
% plot the time evolution of the wind torque
fig_Mr = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector, Mr, 'LineWidth', line_width)
hold on
plot(delta_t_vector, Mg, '--', 'LineWidth', line_width)
hold off
xlabel('Time [s]')
ylabel('Wind and generator torques [Nm]')
title('Wind torque')
legend('Rotor torque', 'Generator torque', 'Location','southeast')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_Mr, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_Mr.png'],'png');

%% Q5: CALCULATION OF THE SLIP
omega_nom = omega(end); % nominal velocity (rad/s)
SL = (omega_nom - sincronous_velocity) / sincronous_velocity; % syncronous motor slip
disp(strcat('The slip is ', num2str(SL), ' = ', num2str(SL*100), ' %'))
disp(['Attention to the fact that omega_nom is estimated from the ' ...
  'results of previous point, and so must be clearly stated '])

omega_item = 30;
V0_q5_vect = linspace(4.5, 20, omega_item);
omega_vect = zeros(omega_item, 1);
P_q5 = zeros(omega_item, 1);

for i=1:omega_item
  V0_q5 = V0_q5_vect(i)*ones(1, i_max);
  % compute the torque and the velocity
  [~, Mg_q5, omega_q5, ~] = control_function(i_max, ...
    omega_initial, R, V0_q5, rho, cP_vs_lambda, MG_coefficient, ...
    sincronous_velocity, delta_t, I);
  omega_vect(i) = omega_q5(end);
  P_q5(i) = omega_vect(i)*Mg_q5(end);
end

fig_omega_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_q5_vect, omega_vect, 'LineWidth',line_width)
grid on
xlabel('V0 [m/s]')
ylabel('\omega [rad/s]')
title('Rotational speed for different velocities')
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_omega_V0, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_omega_V0.png'],'png');

fig_power_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_q5_vect, P_q5/1000, 'LineWidth', line_width)
grid on
xlabel('V0 [m/s]')
ylabel('Generator power P_G [kW]')
title('Generator power for different velocities')
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_power_V0, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_power_V0.png'],'png');

%% Q6: COMPUTATION OF THE ANGULAR VELOCITY, ROTOR AND GENERATOR TORQUE FOR 
%% TURBOLENT WIND INFLOW

clear lambda
V0_t = V0_vs_time(:,2); % load velocity (m/s)
V0_item = size(V0_t, 1); % number of item of the velocity vector
t_max = V0_vs_time(end,1); % maximum time of the recorded data
delta_t_t = t_max / V0_item; % delta_t for the turbolent case (s)

[Mr_t, Mg_t, omega_t, delta_t_vector_t] = control_function(V0_item, ...
  omega_initial, R, V0_t, rho, cP_vs_lambda, MG_coefficient, ...
  sincronous_velocity, delta_t_t, I);

% plot the time evolution of the roataional speed
fig_omega_t = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector_t, omega_t, 'LineWidth', line_width)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
title('Rotor angular velocity - Turbolent wind')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_omega_t, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_omega_t.png'],'png');


% plot the time evolution of the generator torque
fig_Mg_t = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector_t, Mg_t, 'LineWidth', line_width)
xlabel('Time [s]')
ylabel('Generator torque [Nm]')
title('Generator torque - Turbolent wind')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_Mg_t, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_Mg_t.png'],'png');


% plot the time evolution of the wind torque
fig_Mr_t = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector_t, Mr_t, 'LineWidth', line_width)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque - Turbolent wind')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_Mr_t, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_Mr_t.png'],'png');


% plot the time evolution of both torques on the same graphs
fig_Mg_Mt_t = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector_t, Mr_t, 'LineWidth', line_width)
hold on
plot(delta_t_vector_t, Mg_t, '--', 'LineWidth', line_width*0.5)
hold off
legend('Rotor', 'Generator', 'Location','south')
xlabel('Time [s]')
ylabel('Torque [Nm]')
title('Rotor and generator torque')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_Mg_Mt_t, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_Mg_Mt_t.png'],'png');


%% SIMULATION OF THE FREE RUNNING TURBINE

[Mr_f, Mg_f, omega_f, delta_t_vector_f] = control_function(i_max, ...
  omega_initial_f, R, V0_q1, rho, cP_vs_lambda, 0, ...
  sincronous_velocity, delta_t_t, I);

% plot the time evolution of the roataional speed
fig_omega_free = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector_f, omega_f, 'LineWidth', line_width)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
title('Rotor angular velocity - Free running')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_omega_free, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_omega_free.png'],'png');

% plot the time evolution of the wind torque
fig_torque_free = figure('Position', get(0, 'Screensize'));
plot(delta_t_vector_f, Mr_f, 'LineWidth', line_width)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque - Free running')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
saveas(fig_torque_free, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\exercise_control\images\fig_torque_free.png'],'png');

disp(['The velocity reached is the one for which lambda is such that cp=0.' ...
  'This means that the rotor accelerates up to reacing an equilibrium and ' ...
  'then stops'])