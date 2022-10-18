clear 
close all
clc

%%
r_vector = [200 300 400 500 600 645 690]/1000; % (m) radius section
r_size = size(r_vector, 2);
beta_vector = [23 20 16 14.5 13 12.3 11.7]; % (deg) twist angle
chord_vector = [106 117 112 103 88 82 0]/1000; % (m) chord lenght
R = max(r_vector); % blade length
V0_vector = 1:1:40; % velocity vector
V0_item = size(V0_vector, 2);
cd = 0.008;
rho = 1.225; % air density (kg/m^3)
a_guess = 0; % initial guess for a
a_prime_guess = 0; % initial guess for a_prime
fake_zero = 1e-8;
i_max = 1000;
B = 2; % number of blades
omega = pi*2600/30;
beta_correction = 0.1;
T = zeros(V0_item);
P = zeros(V0_item);

for j = 1:V0_item
  dT = zeros(r_size, 1);
  dM = zeros(r_size, 1);
  V0 = V0_vector(j);
  lambda = omega*R/V0;
  % Exercise 1: Implementation of the BEM code
  for i = 1:r_size - 1
    r = r_vector(i);
    beta = beta_vector(i);
    c = chord_vector(i);
    sigma = sigma_function(c, B, r);
    

    cl = 0.1*beta + 0.4;

    [a, a_prime, ct, cn, phi, ~] = induction_factor_convergence(a_guess, ...
      a_prime_guess, R, r, lambda, cd, cl, B, sigma, fake_zero,  i_max, ...
      beta_correction);
    
    dT(i) = B*0.5*rho*c*((1 + a)*V0/sin(phi))^2*cn;
    dM(i) = B*r*0.5*rho*c*(1 + a)*V0*(1 - a_prime)*omega*r*ct/(sin(phi)* ...
      cos(phi));
  
  end
  
  dT(r_size) = 0;
  dM(r_size) = 0;
  [T(j), M] = trapezoidal_integral(dT, dM, r_vector, r_size);
  P(j) = M*omega;
end

figure()
plot(V0_vector, T)
title('Thrust')
xlabel('Wind speed (m/s)')
ylabel('Thrust (N)')

figure()
plot(V0_vector, P/1000)
title('Power (KW)')
xlabel('Wind speed (m/s)')
ylabel('Power (W)')


%% Velocity of the boat
% No acceleration when the trust is equal to the thrust force and so when
