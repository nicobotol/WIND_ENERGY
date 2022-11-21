function [cos_phi] = transformer2( Vb, Vpoc, Cc, Rc, Lc, omega, ...
  L2_prime, R2_prime, n, Lm, R1, L1, Pb )

Vpoc_prime = Vpoc*n ; % Vpoc
Vpoc_prime = Vpoc_prime / sqrt(3); % transform the voltage to phase

Cc_prime = n^2*Cc;
Rc_prime = n^2*Rc;
Lc_prime = n^2*Lc;
% 
% B = [R1+1j*omega*L1 1j*omega*Lm 0 0 0;
%   1 -1 -1 0 0;
%   0 1j*omega*Lm -(R2_prime + Rc_prime + 1j*omega*(L2_prime + Lc_prime)) 0 -1/(1j*omega*Cc_prime);
%   0 0 -1 1 1;
%   0 0 0 0 1/(1j*omega*Cc_prime)];
% 
% 
% res = B \ [Vb 0 0 0 Vpoc_prime]';
% Ib = res(1)
% 
% phi = atan(imag(Ib) / real(Ib));

alpha = (R2_prime + Rc_prime)^2 + omega^2*(L2_prime + Lc_prime)^2;
beta = -omega*(L2_prime + Lc_prime)*Vb;
gamma = (omega*(L2_prime + Lc_prime).*Pb*R1)/Vb;
delta = (omega*L1.*Pb)*(R2_prime + Rc_prime)/Vb;
epsilon = omega*Vpoc_prime*(L2_prime + Lc_prime);
kappa = omega*Cc_prime*Vpoc_prime;

cos_phi = (gamma - delta)/(kappa*alpha - epsilon - beta);

end