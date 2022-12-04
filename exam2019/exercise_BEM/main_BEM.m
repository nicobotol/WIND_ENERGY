% Exercise 1: Implementation of the BEM code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this version of the code the airfoil data are considered independent
% and the angle of attack is provided
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
V0 = 8; % wind speed (m/s)
omega = 2.61; % rotational speed (rad/s)
Theta_p = -3*pi/180; % pitch angle (rad)
beta = 2*pi/180; % twist angle (rad)
c = 1.5; % chord lenght (m)
cl = 0.5; % load coefficient
cd = 0.01; % drag coefficient
r = 24.5; % actual radius (m)
R = 31; % blade length (m)
B = 3; % number of blades
rho = 1.225; % air density (kg/m^3)
a_guess = 0; % initial guess for a
a_prime_guess = 0; % initial guess for a_prime
fake_zero = 1e-8;
i_max = 200;

lambda = omega * R / V0;
sigma = sigma_function(c, B, r); % sigma may be compute independently only once for each section r
% compute the a and a_prime with the iterative method
[a, a_prime, ct, cn, phi, F] = induction_factor_convergence(a_guess, a_prime_guess, ...
  R, r, lambda, beta, Theta_p, cd, cl, c, B, sigma, fake_zero, i_max);

V_rel = sqrt( (omega*r*(1 + a_prime))^2 + (V0*(1 - a))^2);

pn = 0.5 * rho * V_rel^2 * c * cn;
pt = 0.5 * rho * V_rel^2 * c * ct;

display(strcat('a = ', num2str(a)))
display(strcat('a_prime = ', num2str(a_prime)))
display(strcat('pn = ', num2str(pn), '[N/m]'))
display(strcat('pt = ', num2str(pt), '[N/m]'))
display(strcat('F = ', num2str(F)))