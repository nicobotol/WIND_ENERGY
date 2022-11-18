function [Ppoc, Qpoc] = transformer(Ib, Vb, Vpoc, Cc, Rc, Lc, omega, ...
  L2_prime, R2_prime, n, Lm, R1, L1 )

Vpoc_prime = Vpoc*n ; % Vpoc
Vpoc_prime = Vpoc_prime / sqrt(3); % transform the voltage to phase

Cc_prime = n^2*Cc;
Rc_prime = n^2*Rc;
Lc_prime = n^2*Lc;

B = [1 1 0 0;
  1j*omega*Lm 0 0 0;
  1j*omega*Lm -(R2_prime+Rc_prime+1j*omega*(L2_prime+Lc_prime)) -1/(1j*omega*Cc_prime) 0;
  0 -1 1 1];

p_max = size(Ib, 2);
Ipoc = zeros(p_max, 1);
for p=1:p_max
res = B \ [Ib(p) (Vb-Ib(p)*(R1+1j*omega*L1)) 0 0 ]';
Ipoc(p) = res(4);
end

Spoc = 3*Vpoc_prime*conj(Ipoc);
Ppoc = real(Spoc);
Qpoc = imag(Spoc);

end