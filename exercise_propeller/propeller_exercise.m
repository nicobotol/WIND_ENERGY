clear 
close all
clc
%% Load airfoil data
r_vector = [200 300 400 500 600 645 690]/1000; % (m) radius section
r_size = size(r_vector, 2);
beta_vector = [23 20 16 14.5 13 12.3 11.7]; % (deg) twist angle
chord_vector = [106 117 112 103 88 82 0]/1000; % (m) chord lenght
V0_size = 100;
V0_vector = linspace(1, 40, V0_size);
line_width = 1.25;

%% Compute power and torque to drive the propeller

omega = 2600*pi/30; % angular velocity, rad/s
fake_zero = 1e-8; % fake zero to test convergence
beta_correction = 0.1; % correction for compute a and a_prime
i_max = 500; % max number of iterations
rho = 1.225; % air density
B = 2; % number of blades
cd_constant = 0.008;
R = max(r_vector); % propeller radius

T = zeros(V0_size,1);
P = zeros(V0_size,1);

for j=1:V0_size
  V0 = V0_vector(j); % chose the velocity to investigate
  dT = zeros(r_size,1);
  dP = zeros(r_size,1);
  a_guess = 0; % initialize a
  a_prime_guess = 0; % initialize a_prime

  for i = 1:r_size-1
    
    c = chord_vector(i);
    r = r_vector(i); % radius
    theta = beta_vector(i);

    [a, a_prime, phi, cn, ct] = induction_factor(a_guess, a_prime_guess, ...
      omega, r, R, V0, theta, B,cd_constant, c, ...
      beta_correction, fake_zero, i_max );
  
    dT(i) = 1/2*B*rho*c*(1 + a)^2*V0^2/sin(phi)^2*cn; % compute thrust
    dM = 1/2*B*r*rho*c*(1 + a)*V0*(1 - a_prime)*omega*r*ct/(sin(phi)* ...
      cos(phi)); % compute torque
    dP(i) = dM*omega; % compute power

  end
  
  % integrate thrust and power
  [T(j), P(j)] = trapezoidal_integral(dT, dP, r_vector, r_size);
end

% plot the results
figure()
subplot(1, 2, 1)
plot(V0_vector, T, 'LineWidth', line_width);
xlabel('Speed [m/s]')
ylabel('Thrust [N]')
title('Thrust to drive the propeller')
grid on

subplot(1, 2, 2)
plot(V0_vector, P, 'LineWidth', line_width)
xlabel('Speed [m/s]')
ylabel('Power [W]')
title('Power to drive the propeller')
grid on

%% Is the motor strong enough to make the propeller rotate at 2600 rpm?
% Looking at the power courve, it could be seent that it stays always below
% the power of 22 (kW), so the motoe provides enough power to achieve this
% objective
%

%% Compute the velocity of swamp boat
drag = 2.44*V0_vector.^2; % drag for different velocities

% velocity of the swamp boat is reached when the drag is eqaul to the
% thrust

for i = 1:V0_size
  if abs(drag(i) - T(i)) < 10
    velocity = sqrt(T(i) / 2.44);
  end 
end

figure()
plot(V0_vector, T, 'LineWidth', line_width);
hold on
plot(V0_vector, drag, 'LineWidth', line_width);
hold off
xline(velocity, '-', {num2str(velocity)})
xlabel('Speed [m/s]')
ylabel('[N]')
legend('Thrust', 'Drag')


disp(strcat('The velocity of the boat is V = ', num2str(velocity), ' [m/s]'))
