function [Theta_y, Theta_z, u_y, u_z, Ky, Kz] = deflection(V0, r_item)
% This function computes the angle and deflaction of the blade

% load data for the static deflaction
static_blade_data = readmatrix("bladestruc.txt");

% load parameters
parameters

% load the blade data from "bladedat.txt"
[r_vector, c_vector, beta_vector, thick_vector] = load_blade_data(blade_filename);

% initialize variables
Ty = zeros(r_item, 1);
Tz = zeros(r_item, 1);
My = zeros(r_item, 1);
Mz = zeros(r_item, 1);
M1 = zeros(r_item, 1);
M2 = zeros(r_item, 1);
Ky = zeros(r_item, 1);
Kz = zeros(r_item, 1);
Theta_y = zeros(r_item, 1); 
Theta_z = zeros(r_item, 1); 
u_y = zeros(r_item, 1); 
u_z = zeros(r_item, 1); 

% step 0 - internal loads
[py, pz, Theta_p] = internal_beam_loads(V0);

% step 1 - internal forces
for i = r_item:-1:2
  x_diff = r_vector(i) - r_vector(i-1);
  Ty(i-1) = Ty(i) + 0.5*(py(i) + py(i-1))*x_diff;
  Tz(i-1) = Tz(i) + 0.5*(pz(i) + pz(i-1))*x_diff;
  My(i-1) = My(i) - Tz(i)*x_diff - (pz(i-1)/6 + pz(i)/3)*x_diff^2;
  Mz(i-1) = Mz(i) + Ty(i)*x_diff + (py(i-1)/6 + py(i)/3)*x_diff^2;
end

% step 2 - curvature
for i = 1:r_item
  angle = beta_vector(i) + static_blade_data(i, 2)*pi/180 + Theta_p;
  EI1 = static_blade_data(i, 4);
  EI2 = static_blade_data(i, 5);
  M1 = My(i)*cos(angle) - Mz(i)*sin(angle);
  M2 = My(i)*sin(angle) + Mz(i)*cos(angle);
  K1 = M1 / (EI1);
  K2 = M2 / (EI2);
  Ky(i) = K1*cos(angle) + K2*sin(angle);
  Kz(i) = -K1*sin(angle) + K2*cos(angle);
end

% step 3 - angle and deflection
for i = 1:r_item-1
  x_diff = r_vector(i+1) - r_vector(i);
  Theta_y(i+1) = Theta_y(i) + 0.5*(Ky(i+1) + Ky(i))*x_diff;
  Theta_z(i+1) = Theta_z(i) + 0.5*(Kz(i+1) + Kz(i))*x_diff;
  u_y(i+1) = u_y(i) + Theta_z(i)*x_diff + (Kz(i+1)/6 + Kz(i)/3)*x_diff^2;
  u_z(i+1) = u_z(i) - Theta_y(i)*x_diff - (Ky(i+1)/6 + Ky(i)/3)*x_diff^2;
end
end