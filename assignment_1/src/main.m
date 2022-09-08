clear all 
close all
clc

parameters
filename = 'airfoil_data\FFA-W3-241.txt';
FFAW3241 = readmatrix(filename)

% lambda = omega * R / V0;
% % compute the a and a_prime with the iterative method
% [a, a_prime, ct, cn, phi] = induction_factor_convergence(a, a_prime, R, r, lambda, beta, Theta_p, cd, cl, c, B, fake_zero, i_max);
% 
% a 
% a_prime
% 
% V_rel = sqrt( (omega*r*(1 + a_prime))^2 + (V0*(1 - a))^2);
% 
% pn = 0.5 * rho * V_rel^2 * c * cn 
% pt = 0.5 * rho * V_rel^2 * c * ct

