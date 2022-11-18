function [Ig, Vg, fg, omega_E, delta] = rotor_equation(Rs, Ls, pp, flux, omega_M, Pmech)
% This function solves the rotor electrical equation for each mechanical
% rotational speed

omega_E = pp*omega_M; % electrical rotational speed
fg = omega_E/(2*pi); % electrical frequency (Hz)

delta = 0.5.*asin(2/3*Pmech*Ls ./ (flux.^2.*omega_E)); % angle between the back 
% emf and Vg
Ig = Pmech./(flux.*omega_E.*cos(delta)*3); % generator current (A)
Vg = flux.*omega_E.*cos(delta) - Ig*Rs; % generator tension (V)

end