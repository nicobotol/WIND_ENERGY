clear all
close all
clc

%% CONSTANT PARAMETERS
time_max = 50; % maximum time to be studied (s)
delta_t = 0.05; % discretization time (s)
i_max = time_max / delta_t; % number of iterations

sincronous_velocity = 5; % (rad/s)
MG_coefficient = 2.55e5; % generator curve coefficient

R = 9.5; % rotor diameter (m)
rho = 1.225; % air density (kg/m^3)
V0 = 10; % wind speed (m/s)
omega_initial = 1; % initial angular velocity (rad/s)
delta_t = 0.01; % time step (s)
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
%% COMPUTATION OF THE ANGULAR VELOCITY
Mr = zeros(i_max, 1);
Mg = zeros(i_max, 1);
omega = zeros(i_max, 1);
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
plot(delta_t_vector, Mg)
xlabel('Time [s]')
ylabel('Wind torque [Nm]')
title('Wind torque')
