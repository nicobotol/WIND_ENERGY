clear all 
close all
clc

% load the physical parameters from the file "parameters.m"
parameters

% load the airfoil data from the different files
[aoa, cl_mat, cd_mat] = load_airfoil_data(filenames);

% load the blade data from "bladedat.txt"
[r, c, beta, thick] = load_blade_data(blade_filename);

% function to interpolate cl and cd for the given value of alpha and t/c
[cl, cd] = interpolate_cl_cd(aoa, cl_mat, cd_mat, thick_prof, -130, thick(18));

lambda = omega * R / V0;
% compute the a and a_prime with the iterative method
[a, a_prime, ct, cn, phi] = induction_factor_convergence(a, a_prime, R, r, lambda, beta, Theta_p, cd, cl, c, B, fake_zero, i_max);

V_rel = sqrt( (omega*r*(1 + a_prime))^2 + (V0*(1 - a))^2);

pn = 0.5 * rho * V_rel^2 * c * cn;
pt = 0.5 * rho * V_rel^2 * c * ct;
