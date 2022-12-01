function [Ppoc, Qpoc] = Q3_new(P_g, Z1, Z2_prime, Zm, Zc, Vpoc_prime,  Z_cable)
  
Pa = P_g / 3; % phase power on the generator side [W]
Pb = Pa;
Qb = 0.2*Pb;
% phi = atan(Qb/Pb); % leading angle
% Sb_mag = sqrt(Pb.^2 + Qb.^2);
Sb = Pb + 1j*Qb; % complex power on b side (VAR)

syms I2 Im

eq1 = Vpoc_prime - (Z2_prime + Z_cable)*I2 - Im*Zm == 0;
eq2 = Sb/(conj(Im - I2) ) - Z1*(-I2 + Im) - Im*Zm == 0;

sol = vpasolve([eq1, eq2], [I2, Im]);

Ic = Vpoc_prime/Zc;

Ipoc = subs(Ic + I2, sol);
Spoc = 3*Vpoc_prime*conj(Ipoc);

Ppoc = real(Spoc);
Qpoc = imag(Spoc);

end