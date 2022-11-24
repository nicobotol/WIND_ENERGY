clc
% clear all

parameters;

Pa = 322.233516924289;
Pb = Pa/sqrt(3);
Qb = 0.2*Pb;
Sb = Pb + 1j*Qb;
omega = 314;

Cc_prime = n^2*Cc;
Rc_prime = n^2*Rc;
Lc_prime = n^2*Lc;
  Z1 = R1 + 1j*omega*L1;
Z2_prime = (R2_prime + Rc_prime) + 1j*omega*(L2_prime + Lc_prime);
Zm =  1j*omega*Lm;
Zc = 1/(1j*omega*Cc_prime);

for p=1:47
Pa = P_g(p);
Pb = Pa/sqrt(3);
Qb = 0.2*Pb;
Sb = Pb + 1j*Qb;
for i = 1:100
  Vb = 4000/sqrt(3)*exp(1j*pi/180*i/100);
  Ib = conj(Sb/Vb);

  Vab = Vb - Ib*Z1;
  Im = Vab/Zm;
  I2 = Ib - Im;
  Vpoc_prime = Vab - I2*Z2_prime;
  disp(abs(Vpoc_prime))
  if abs(Vpoc_prime) - 33000/sqrt(3)*n < 0.001 
    disp('convergence ok')
    break
  else
  end 

end
Ic = Vpoc_prime ./ Zc;
Ipoc = I2 - Ic;

Spoc = 3*Vpoc_prime.*conj(Ipoc); % the sign (-) comes from the fact
Ppoc(p) = real(Spoc)
Qpoc(p) = imag(Spoc)
end