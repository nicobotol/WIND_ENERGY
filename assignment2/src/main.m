%% Wind Turbine Technology and Aerodynamics - Assignment 2

clear 
close all
clc

%% Load data of the problem
parameters

%% Question 1
omega_M = 0.001:0.5:omega_M_max; % (rpm)
omega_M = omega_M*pi/30; % (rad/s)

Pmech = 0.5*rho*A*cp_opt*(omega_M*R/lambda_opt).^3; % mechanical power from the wind (W)
[Ig, Vg, fg, omega_E, delta] = rotor_equation(Rs, Ls, pp, flux, omega_M, Pmech);

figure()
plot(omega_M, Pmech)
xlabel('Rotational speed [rad/s]')
ylabel('Incoming power [W]')

figure()
plot(omega_M, Ig)
xlabel('Rotational speed [rad/s]')
ylabel('Ig [A]')

figure()
plot(omega_M, Vg)
xlabel('Rotational speed [rad/s]')
ylabel('Vg [V]')

figure()
plot(omega_M, fg)
xlabel('Rotational speed [rad/s]')
ylabel('Generator frequency [Hz]')

%% Question 2

% Losses
P_loss = Ig.^2*Rs; % Loss in the resistance (W)
Q_loss = Ig.^2.*omega_E*Ls; % Loss in the inductance (VAR)

% Output power
P_g = flux.*omega_M.*cos(delta); % (W)
Q_g = flux.*omega_M.*sin(delta); % (VAR)
S_g = sqrt(P_g.^2 + Q_g.^2); % (VA)

figure()
plot(omega_M, P_loss)
hold on
plot(omega_M, Q_loss)
hold off
legend('P', 'Q')
xlabel('Rotational speed [rad/s]')
ylabel('Power losses')
