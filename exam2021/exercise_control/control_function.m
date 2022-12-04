function [Mr, Mg, omega, delta_t_vector] = control_function(i_max, ...
  omega_initial, R, V0, rho, cP_vs_lambda, MG_coefficient, ...
  sincronous_velocity, delta_t, I)
%UNTITLED3 Summary of this function goes here
% i_max max number of iterations
% omega_initial initial rotational speed
% V0 vector of wind velocity

Mr = zeros(i_max, 1); % vector for rotor torque (Nm)
Mg = zeros(i_max, 1); % vector for generator torque (Nm)
omega = zeros(i_max, 1); % vector for angular velocity (rad/s)
delta_t_vector = zeros(i_max, 1);

% initialize the first elements of the vectors of omega, wind torque and
% generator torque
omega(1,1) = omega_initial;
lambda = omega(1,1) * R / V0(1); % compute lambda
cP = interp1(cP_vs_lambda(:,1), cP_vs_lambda(:,2), lambda); % interpolate the table
Mr(1,1) = MR(rho, V0(1), R, cP, omega(1,1)); % wind torque (Nm)
Mg(1,1) = MG(omega(1,1), MG_coefficient, sincronous_velocity); % generator torque (Nm)

for i=2:i_max
  lambda = omega(i - 1, 1) * R / V0(i - 1); % compute lambda
  cP = interp1(cP_vs_lambda(:,1), cP_vs_lambda(:,2), lambda); % interpolate the table
  %omega_old = omega(i - 1, 1); % store previous omega
  Mr(i,1) = MR(rho, V0(i - 1), R, cP, omega(i - 1, 1)); % wind torque (Nm)
  Mg(i,1) = MG(omega(i - 1, 1), MG_coefficient, sincronous_velocity); % generator torque (Nm)
  %omega(i,1) = omega(i - 1, 1) + delta_t/I*(Mr(i - 1 ,1) - Mg(i - 1,1)); % update omega
  omega(i,1) = omega(i - 1, 1) + delta_t/I*(Mr(i ,1) - Mg(i,1)); % update omega

  delta_t_vector(i, 1) = delta_t_vector(i - 1, 1) + delta_t; % store the time for the plot
end


end