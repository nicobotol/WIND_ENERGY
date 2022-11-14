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
[P_loss, Q_loss, S_loss, P_g, Q_g, S_g] = losses(Ig, omega_E, delta, Ls, Rs, flux);

figure()
plot(omega_M, P_loss)
hold on
plot(omega_M, Q_loss)
hold off
legend('P', 'Q')
xlabel('Rotational speed [rad/s]')
ylabel('Power losses')

% Output complex power
figure()
plot(omega_M, S_g)
xlabel('Rotational speed [rad/s]')
ylabel('Complex power output [VAR]')
%%
% verify wheter the sum of the power is consistent
figure()
plot(omega_M, S_loss + S_g)
hold on
plot(omega_M, Pmech)
hold off
legend('Ploss+Pg', 'Pmech')
xlabel('Rotational speed [rad/s]')
ylabel('Power [W]')

%% Question 3

omega = 2*pi*f_grid; % angular velocity grid side
Vb = Vb_ll / sqrt(3); % phase voltage (V)

[Ppoc, Qpoc] = transformer(P_g, Vb, Cc, n, R1, Lm, L1, omega, Vpoc);

% POC active power
figure()
plot(omega_M, Ppoc)
xlabel('Rotational speed [rad/s]')
ylabel('POC power [W]')