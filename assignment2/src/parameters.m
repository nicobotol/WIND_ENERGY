%% Data of the problem
rho = 1.225; % air density (kg/m^3)
Rs = 0.064; % wtg internal resistnace (ohm)
Ls = 0.0018; % wtg internal inductance (H)
flux = 19.49; % rototr magnetic flux (Wb)
poles = 320; % number of poles of the rotor
pp = poles/2; % pole pairs of the rotor
omega_M_max = 30*1.01/pi; % maximum rotational speed (rpm)
R = 89.17; % generator rotor radius (m)
A = pi*R^2; % blade swapt area
cp_opt = 0.465; % optimal cp
lambda_opt = 7.857; % optimal lambda