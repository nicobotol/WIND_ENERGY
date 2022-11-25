function [Vb_solution, Vab_solution, Vcd_solution, Vpoc_solution, Ppoc, Qpoc] = Q3(P_g, Z1, Z2_prime, Zm, Zc, Vpoc_prime,  Zcable)
%UNTITLED Summary of this function goes here

Pa = P_g / 3; % phase power on the generator side [W]
Pb = Pa;
Qb = 0.2*Pb;
% phi = atan(Qb/Pb); % leading angle
% Sb_mag = sqrt(Pb.^2 + Qb.^2);
Sb = Pb + 1j*Qb; % complex power on b side (VAR)

Vb_vector = linspace(2200, 2500, 1000);
threshold = 5;

for i = 1:1000
  Vb = Vb_vector(i);  
  Ib = conj(Sb / Vb);
  Vab = Vb - Ib*Z1;
  Im = Vab/Zm;
  I2 = Ib - Im;
  Vcd = Vab - Z2_prime*I2;
  Vpoc_star = Vcd - I2*Zcable;
  Vpoc_star_mag = abs(Vpoc_star);

  if( abs(Vpoc_star_mag - Vpoc_prime) <  threshold)
    
    disp(strcat('convergence at', num2str(i)))
    Vb_solution = Vb; 
    Vab_solution = Vab;
    Vcd_solution = Vcd;
    Vpoc_solution = abs(Vpoc_star);
    Ic = Vpoc_star / Zc;
    Ipoc = I2 - Ic;
    Spoc = 3*Vpoc_star*conj(Ipoc);
    %Ppoc = imag(Spoc);
    Ppoc = real(Spoc);
    Qpoc = imag(Spoc);
%     Qpoc = real(3*(Vcd)*conj(I2));
    break
  end

end

end