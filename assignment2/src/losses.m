function [P_loss, Q_loss, S_loss, P_g, Q_g, S_g] = losses(Ig, omega_E, delta, Ls, Rs, flux)
%% This function compute the losses

P_loss = Ig.^2*Rs; % Loss in the resistance (W)
Q_loss = Ig.^2.*omega_E*Ls; % Loss in the inductance (VAR)
S_loss = sqrt(P_loss.^2 + Q_loss.^2);

% Output power
P_g = 3*flux.*omega_E.*Ig.*cos(delta); % (W)
Q_g = 3*flux.*omega_E.*Ig.*sin(delta); % (VAR)
S_g = sqrt(P_g.^2 + Q_g.^2); % (VA)

end