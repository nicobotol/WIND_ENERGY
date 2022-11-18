function [phi] = transformer2( Vb, Vpoc, Cc, Rc, Lc, omega, ...
  L2_prime, R2_prime, n, Lm, R1, L1 )

Vpoc_prime = Vpoc*n ; % Vpoc
Vpoc_prime = Vpoc_prime / sqrt(3); % transform the voltage to phase

Cc_prime = n^2*Cc;
Rc_prime = n^2*Rc;
Lc_prime = n^2*Lc;

B = [1 0 0 0 0;
  0 1 -1 0 0 ;
  -1/(1j*omega*Cc_prime) -(R2_prime+Rc_prime+1j*omega*(L2_prime+Lc_prime)) 0 1j*omega*Lm 0;
  0 0 0 1j*omega*Lm R1+1j*omega*L1;
  0 1 0 1 1];


res = B \ [1j*Vpoc_prime*omega*Cc_prime 1j*Vpoc_prime*omega*Cc_prime 0 Vb 0]';
Ib = res(5)

phi = atan(imag(Ib) / real(Ib));

end