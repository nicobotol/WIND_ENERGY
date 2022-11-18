function [P_loss, Q_loss, S_loss, P_g, Q_g, S_g] = losses(Ig, Vg, omega_E, Ls, Rs)
%% This function compute the losses

P_loss = 3*Ig.^2*Rs; % Loss in the resistance (W)
Q_loss = 3*Ig.^2.*omega_E*Ls; % Loss in the inductance (VAR)
S_loss = sqrt(P_loss.^2 + Q_loss.^2);

% Output power
S_g = 3*Vg.*conj(Ig); % complex output power
P_g = real(S_g);
Q_g = imag(S_g);

end