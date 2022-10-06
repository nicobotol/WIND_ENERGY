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

%% Q1 & Q2: PLOT OF THE LOADED DATA
% load the data form the table
cP_vs_lambda = load('Cp_Lambda.txt');
V0_vs_time = load('usim.dat');

figure()
plot(cP_vs_lambda(:,1), cP_vs_lambda(:,2))
xlabel('\lambda')
ylabel('cP')
title('cP as function of \lambda')
grid on

figure()
plot(V0_vs_time(:,1), V0_vs_time(:,2))
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
figure()
plot(delta_t_vector, omega)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
title('Rotor angular velocity')
grid on

% plot the time evolution of the generator torque
figure()
plot(delta_t_vector, Mg)
xlabel('Time [s]')
ylabel('Generator torque [Nm]')
title('Generator torque')
grid on

% plot the time evolution of the wind torque
figure()
plot(delta_t_vector, Mr)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque')
grid on

%% Q5: CALCULATION OF THE SLIP
omega_nom = omega(end); % nominal velocity (rad/s)
SL = (omega_nom - sincronous_velocity) / sincronous_velocity; % syncronous motor slip

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
figure()
plot(delta_t_vector_t, omega_t)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
title('Rotor angular velocity - Turbolent wind')
grid on

% plot the time evolution of the generator torque
figure()
plot(delta_t_vector_t, Mg_t)
xlabel('Time [s]')
ylabel('Generator torque [Nm]')
title('Generator torque - Turbolent wind')
grid on

% plot the time evolution of the wind torque
figure()
plot(delta_t_vector_t, Mr_t)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque - Turbolent wind')
grid on

% plot the time evolution of both torques on the same graphs
figure()
plot(delta_t_vector_t, Mr_t)
hold on
plot(delta_t_vector_t, Mg_t)
hold off
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque')
grid on

%% SIMULATION OF THE FREE RUNNING TURBINE

[Mr_f, Mg_f, omega_f, delta_t_vector_f] = control_function(i_max, ...
  omega_initial_f, R, V0_q1, rho, cP_vs_lambda, 0, ...
  sincronous_velocity, delta_t_t, I);
% plot the time evolution of the roataional speed
figure()
plot(delta_t_vector_f, omega_f)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
title('Rotor angular velocity - Free running')
grid on

% plot the time evolution of the wind torque
figure()
plot(delta_t_vector_f, Mr_f)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque - Free running')
grid on