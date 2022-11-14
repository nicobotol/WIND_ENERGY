function [Ig, Vg, fg, omega_E, delta] = rotor_equation(Rs, Ls, pp, flux, omega_M, Pmech)
% This function solves the rotor electrical equation for each mechanical
% rotational speed

omega_E = pp*omega_M; % electrical rotational speed
fg = omega_E/(2*pi); % electrical frequency (Hz)

delta = 0.5.*asin(2.*Pmech.*omega_E*Ls ./ (flux.*omega_M).^2); % angle between the back 
% emf and Vg
Ig = sqrt(3).*Pmech./(flux.*omega_M.*cos(delta)*3); % generator current (A)
Vg = flux.*omega_M.*cos(delta)/sqrt(3) - Ig*Rs; % generator tension (V)

end