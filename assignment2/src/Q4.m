function [delta] = Q4(P_g, Z1, Z2_prime, Zm, Zc, Vpoc_prime,  Zcable)
%UNTITLED Summary of this function goes here

Pa = P_g / 3; % phase power on the generator side [W]

Ipoc_vector = linspace(-10000, 10000, 10000);
threshold = 1000;

for i = 1:10000
  Ipoc = Ipoc_vector(i);
  Ic = Vpoc_prime / Zc;
  I2 = Ipoc - Ic;
  Vab = Vpoc_prime - I2*(Z2_prime + Zcable);
  Im = Vab/Zm;
  Ib = I2 - Im;
  Vb = Vab - Ib*Z1;
  Pg = abs(real(Vb*conj(Ib)));

  if( abs(Pg - Pa) >  threshold)
    disp(strcat('convergence at', num2str(i)))
    Qg = imag(Vb*conj(Ib));
    delta = atan(Qg/Pg);
    break
  end

end

end