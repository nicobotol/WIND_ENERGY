function [Ppoc, Qpoc] = transformer(Ib, Vb, Vpoc, Cc, Rc, Lc, omega, L2_prime, R2_prime, n, Lm, R1, L1 )

Vpoc_prime = Vpoc / n; % Vpoc
Vpoc_prime = Vpoc_prime / sqrt(3); % transform the voltage to phase

Cc_prime = n^2*Cc;
Rc_prime = n^2*Rc;
Lc_prime = n^2*Lc;

% Zm = 1j*Lm*omega;
% Zb = R1 + 1j*L1;
% Zc = R2_prime + Rc_prime + 1j*(L2_prime + Lc_prime);
% 
% Ipoc = -(Vb - Ib*Zb)/Zm - Vpoc/Zc + Ib;

p_max = size(Ib, 2);
Ipoc = zeros(p_max, 1);
for p=1:p_max
  res = [1j*omega*Lm R1+1j*omega*L1 0 0 ;
    Lm -(R2_prime + Rc_prime +1j*omega*(L2_prime + Lc_prime)) 0 1j/omega*Cc_prime;
    0 0 0 1j*omega*Cc_prime;
    0 1 -1 0] / [Vb 0 Vpoc_prime Ib(p)];

  Ipoc(p) = res(3);
end

Spoc = Vpoc_prime*conj(Ipoc);
Ppoc = real(Spoc);
Qpoc = imag(Spoc);



% % Vb line-line voltage VSC-B (V)
% 
% Vb = Vb / sqrt(3);
% Cc1 = Cc*n^2; % report the wire capacitor to the primary
% Vpoc1 = Vpoc / (sqrt(3)*n); % report the POC voltage to the primary
% 
% 
% Pb = Pa; % incoming active power
% Qb = 0.2.*Pb; % incoming reactive power
% phi = atan(Qb./Pb); % phase angle incoming quantities
% 
% Ib = sqrt(Pb.^2 + Qb.^2).*exp(-1i.*phi) / Vb;
% 
% Vc = Vb - Ib*(R1 + 1i*omega*L1);
% 
% I3 = Ib.*Vc/(1i*omega*Lm);
% I4 = Vpoc*omega*Cc1/(-1i);
% 
% Ipoc = I3 - I4;
% Spoc = Vpoc1.*conj(Ipoc); % complex power delivered to the POC
% Ppoc = real(Spoc); % real power delivered to the POC
% Qpoc = imag(Spoc); % complex power deliverd to the POC

end