function [Ppoc, Qpoc] = transformer(Pa, Vb, Cc, n, R1, Lm, L1, omega, Vpoc)
% Vb line-line voltage VSC-B (V)

Vb = Vb / sqrt(3);
Cc1 = Cc*n^2; % report the wire capacitor to the primary
Vpoc1 = Vpoc / (sqrt(3)*n); % report the POC voltage to the primary


Pb = Pa; % incoming active power
Qb = 0.2.*Pb; % incoming reactive power
phi = atan(Qb./Pb); % phase angle incoming quantities

Ib = sqrt(Pb.^2 + Qb.^2).*exp(-1i.*phi) / Vb;

Vc = Vb - Ib*(R1 + 1i*omega*L1);

I3 = Ib.*Vc/(1i*omega*Lm);
I4 = Vpoc*omega*Cc1/(-1i);

Ipoc = I3 - I4;
Spoc = Vpoc1.*conj(Ipoc); % complex power delivered to the POC
Ppoc = real(Spoc); % real power delivered to the POC
Qpoc = imag(Spoc); % complex power deliverd to the POC

end