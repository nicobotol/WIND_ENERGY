clear 
close all
clc

omega = 115*pi/30; % [rad/s]
V0 = 12; % wind velocity [m/s]
R = 2; % radius [m]
rho = 1.225; % air density [kg/m^3]
theta = 1.5*pi; % azimuthal angle [rad]
Wx = 4; % induced wind on x direction [m/s]
Wy = 0; % induced wind on y direction [m/s]

[Vrel, alpha] = compute_p_exam2020(V0, R, omega, theta, Wx, Wy);
disp(strcat('The relative velocity is ', num2str(Vrel), ' [m/s]'))
disp(strcat('The angle of attack is ', num2str(alpha), ' [rad]'))