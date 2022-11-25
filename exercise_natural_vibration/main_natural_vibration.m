clear 
close all
clc

%% Load data
% load the physical parameters from the file "parameters.m"
parameters

% load pitchig angle
load('Pitching.mat'); % first row: windspeed between 4-25 (m/s)
                      % second row: angle for stall
                      % third row: angle for feathering

% load the airfoil data from the different files
[aoa_mat, cl_mat, cd_mat] = load_airfoil_data(filenames);

% load the blade data from "bladedat.txt"
[r_vector, c_vector, beta_vector, thick_vector] = load_blade_data(blade_filename);

r_item = size(r_vector, 2); % number of cross sections along the blade
r_item_no_tip = r_item - 1; % numer of cross section to take into account

% build the mass matrix
[M] = build_mass(r_item);

% build the flexibility matrix
F = zeros(2*r_item_no_tip);
for n=3:2*r_item
 [u] = deflection_eigen(r_item, n);  
 F(:, n - 2) = u;
 
end

% solve the eigenvalue problem
[V, D] = eigs(F*M);

% find the eigenvalues
omega1 = sqrt(1/D(1,1))
omega2 = sqrt(1/D(2,2))
omega3 = sqrt(1/D(3,3))
omega4 = sqrt(1/D(4,4))
omega5 = sqrt(1/D(5,5))
omega6 = sqrt(1/D(6,6))

% find the eigenvectors
MS1_y = V(1:2:end, 1)/max(abs(V(:, 1)));
MS1_z = V(2:2:end, 1)/max(abs(V(:, 1)));

MS2_y = V(1:2:end, 3)/max(abs(V(:, 3)));
MS2_z = V(2:2:end, 3)/max(abs(V(:, 3)));

MS3_y = V(1:2:end, 5)/max(abs(V(:, 5)));
MS3_z = V(2:2:end, 5)/max(abs(V(:, 5)));

figure()
plot(r_vector(1:17), MS1_y )
hold on
plot(r_vector(1:17), MS2_y )
plot(r_vector(1:17), MS3_y )