% This code answers to question 4 of the exam 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A wind turbine has a radius of R=30 m and a pitch angle of θp=1 degrees. 
% The free wind speed is Vo=9 m/s, the wind speed at the rotor plane is 
% u=6m/s and the tip speed ratio is λ=7.
% At r=15 m/s the angle of attack is measured to α=2 degrees and the twist
% β=7 degrees.
% Calculate the tangential induction factor a'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close 
clc

% parameters
V0 = 9; % wind speed (m/s)
lambda = 7; % tip speed ratio
R = 30; % blade length (m)
Theta_p = 1*pi/180; % pitch angle (rad)
u = 6; % wind speed at rotor plane (m/s)
beta = 7*pi/180; % twist angle (rad) at distance r
r = 15; % actual radius (m)
aoa = 2*pi/180; % angle of attach at the current balde position

omega = V0*lambda/R; % rotational speed (rad/s)
a = 1 - u/V0; % axial indction factor
phi = aoa + beta + Theta_p; % flow angle (rad)

a_prime = ((1 - a)*V0)/(omega*r*tan(phi)) - 1;

display(strcat('a_prime = ', num2str(a_prime)))
