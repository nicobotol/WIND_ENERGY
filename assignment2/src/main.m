%% Wind Turbine Technology and Aerodynamics - Assignment 2

clear 
close all
clc

%% Load data of the problem
parameters

%% Question 1
omega_M = 0.001:0.5:omega_M_max; % (rpm)
omega_M = omega_M*pi/30; % (rad/s)
rpm = 30/pi*omega_M;

Pmech = 0.5*rho*A*cp_opt*(omega_M*R/lambda_opt).^3; % mechanical power from the wind (W)
[Ig, Vg, fg, omega_E, delta] = rotor_equation(Rs, Ls, pp, flux, omega_M, Pmech);

figure()
plot(rpm, Pmech)
xlabel('Rotational speed [rpm]')
ylabel('Incoming power [W]')

figure()
plot(rpm, Ig)
xlabel('Rotational speed [rpm]')
ylabel('Ig [A]')

figure()
plot(rpm, Vg)
xlabel('Rotational speed [rpm]')
ylabel('Vg [V]')

figure()
plot(rpm, fg)
xlabel('Rotational speed [rpm]')
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
clc
Vb = Vb_ll / sqrt(3); % phase voltage on VSC-B (V)

Pa = Vg.*Ig;
Pb = Pa;
Qb = 0.2*Pb;
phi = atan(Qb/Pb); % leading angle
Sb_mag = sqrt(Pb.^2 + Qb.^2);
Sb = Sb_mag*exp(1j*phi); % complex power on b side (VAR)
Ib = conj(Sb / Vb);


omega = 2*pi*f_grid; % angular velocity grid side


[Ppoc, Qpoc] = transformer(Ib, Vb, Vpoc, Cc, Rc, Lc, omega, L2_prime, R2_prime, n, Lm, R1, L1 );

% POC active power
figure()
plot(omega_M, Ppoc)
xlabel('Rotational speed [rad/s]')
ylabel('POC power [W]')