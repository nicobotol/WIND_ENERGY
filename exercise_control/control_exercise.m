clear all
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
V0 = 10; % wind speed (m/s)
V0_f = 10; % wind speed for free running turbine (m/s)
omega_initial = 1; % initial angular velocity (rad/s)
omega_initial_f = 1; % initial angular velocity free running turbine (rad/s)
I = 1.5e4; % rotor inertia (kgm^2)

%% PLOT OF THE LOADED DATA
% load the data form the table
cP_vs_lambda = load('Cp_Lambda.txt');
V0_vs_time = load('usim.dat');

figure()
plot(cP_vs_lambda(:,1), cP_vs_lambda(:,2))
xlabel('\lambda')
ylabel('cP')
title('cP as function of \lambda')

figure()
plot(V0_vs_time(:,1), V0_vs_time(:,2))
title('Turbolent wind')
xlabel('Time (s)')
ylabel('Wind velocity (m/s)')

%% COMPUTATION OF THE ANGULAR VELOCITY, ROTOR AND GENERATOR TORQUE 
Mr = zeros(i_max, 1); % vector for rotor torque (Nm)
Mg = zeros(i_max, 1); % vector for generator torque (Nm)
omega = zeros(i_max, 1); % vector for angular velocity (rad/s)
delta_t_vector = zeros(i_max, 1);

% initialize the first elements of the vectors of omega, wind torque and
% generator torque
omega(1,1) = omega_initial;
lambda = omega(1,1) * R / V0; % compute lambda
cP = interp1(cP_vs_lambda(:,1), cP_vs_lambda(:,2), lambda); % interpolate the table
Mr(1,1) = MR(rho, V0, R, cP, omega(1,1)); % wind torque (Nm)
Mg(1,1) = MG(omega(1,1), MG_coefficient, sincronous_velocity); % generator torque (Nm)

for i=2:i_max
  lambda = omega(i - 1, 1) * R / V0; % compute lambda
  cP = interp1(cP_vs_lambda(:,1), cP_vs_lambda(:,2), lambda); % interpolate the table
  omega_old = omega(i - 1, 1); % store previous omega
  Mr(i,1) = MR(rho, V0, R, cP, omega(i - 1, 1)); % wind torque (Nm)
  Mg(i,1) = MG(omega(i - 1, 1), MG_coefficient, sincronous_velocity); % generator torque (Nm)
  omega(i,1) = omega(i - 1, 1) + delta_t/I*(Mr(i - 1 ,1) - Mg(i - 1,1)); % update omega

  delta_t_vector(i, 1) = delta_t_vector(i - 1, 1) + delta_t; % store the time for the plot
end

% plot the time evolution of the roataional speed
figure()
plot(delta_t_vector, omega)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
title('Rotor angular velocity')

% plot the time evolution of the generator torque
figure()
plot(delta_t_vector, Mg)
xlabel('Time [s]')
ylabel('Generator torque [Nm]')
title('Generator torque')

% plot the time evolution of the wind torque
figure()
plot(delta_t_vector, Mr)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque')

%% CALCULATION OF THE SLIP
omega_nom = omega(end); % nominal velocity (rad/s)
SL = (omega_nom - sincronous_velocity) / sincronous_velocity; % syncronous motor slip

%% COMPUTATION OF THE ANGULAR VELOCITY, ROTOR AND GENERATOR TORQUE FOR 
%% TURBOLENT WIND INFLOW

V0_t = V0_vs_time(:,2); % load velocity (m/s)
V0_item = size(V0_t, 1); % number of item of the velocity vector
t_max = max(V0_vs_time(:,1)); % maximum time of the recorded data
delta_t_t = t_max / V0_item; % delta_t for the turbolent case (s)
Mr_t = zeros(V0_item, 1);
Mg_t = zeros(V0_item, 1);
omega_t = zeros(V0_item, 1);
delta_t_vector_t = zeros(V0_item, 1);

% initialize the first elements of the vectors of omega, wind torque and
% generator torque
omega_t(1,1) = omega_initial;
lambda = omega_t(1,1) * R / V0_t(1); % compute lambda
cP = interp1(cP_vs_lambda(:,1), cP_vs_lambda(:,2), lambda); % interpolate the table
Mr_t(1,1) = MR(rho, V0_t(1), R, cP, omega_t(1,1)); % wind torque (Nm)
Mg_t(1,1) = MG(omega_t(1,1), MG_coefficient, sincronous_velocity); % generator torque (Nm)

for j=2:V0_item
  lambda = omega_t(j - 1, 1) * R / V0_t(j); % compute lambda
  cP = interp1(cP_vs_lambda(:,1), cP_vs_lambda(:,2), lambda); % interpolate the table
  Mr_t(j,1) = MR(rho, V0_t(j), R, cP, omega_t(j - 1, 1)); % wind torque (Nm)
  Mg_t(j,1) = MG(omega(j - 1, 1), MG_coefficient, sincronous_velocity); % generator torque (Nm)
  omega_t(j,1) = omega_t(j - 1, 1) + delta_t_t/I*(Mr(j - 1 ,1) - Mg(j - 1,1)); % update omega

  delta_t_vector_t(j, 1) = delta_t_vector_t(j - 1, 1) + delta_t_t; % store the time for the plot
end

% plot the time evolution of the roataional speed
figure()
plot(delta_t_vector_t, omega_t)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
title('Rotor angular velocity - Turbolent wind')

% plot the time evolution of the generator torque
figure()
plot(delta_t_vector_t, Mg_t)
xlabel('Time [s]')
ylabel('Generator torque [Nm]')
title('Generator torque - Turbolent wind')

% plot the time evolution of the wind torque
figure()
plot(delta_t_vector_t, Mr_t)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque - Turbolent wind')

%% SIMULATION OF THE FREE RUNNING TURBINE

t_max_f = max(V0_vs_time(:,1)); % maximum time of the recorded data
delta_t_f = t_max / V0_item; % delta_t for the turbolent case (s)
Mr_f = zeros(V0_item, 1);
Mg_f = zeros(V0_item, 1);
omega_f = zeros(V0_item, 1);
delta_t_vector_f = zeros(V0_item, 1);

% initialize the first elements of the vectors of omega, wind torque and
% generator torque
omega_f(1,1) = omega_initial_f;
lambda = omega_f(1,1) * R / V0_f; % compute lambda
cP = interp1(cP_vs_lambda(:,1), cP_vs_lambda(:,2), lambda); % interpolate the table
Mr_f(1,1) = MR(rho, V0_f, R, cP, omega_t(1,1)); % wind torque (Nm)

for j=2:V0_item
  lambda = omega_f(j - 1, 1) * R / V0_f; % compute lambda
  cP = interp1(cP_vs_lambda(:,1), cP_vs_lambda(:,2), lambda); % interpolate the table
  Mr_f(j, 1) = MR(rho, V0_f, R, cP, omega_f(j - 1, 1)); % wind torque (Nm)
  omega_f(j, 1) = omega_f(j - 1, 1) + delta_t_f/I*Mr(j - 1 ,1); % update omega

  delta_t_vector_f(j, 1) = delta_t_vector_f(j - 1, 1) + delta_t_f; % store the time for the plot
end

% plot the time evolution of the roataional speed
figure()
plot(delta_t_vector_f, omega_f)
xlabel('Time [s]')
ylabel('\omega [rad/s]')
title('Rotor angular velocity - Free running')

% plot the time evolution of the wind torque
figure()
plot(delta_t_vector_f, Mr_f)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque - Free running')