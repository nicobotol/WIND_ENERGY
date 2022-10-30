function [a_new, a_prime_new, phi, cn, ct] = induction_factor(a_guess, a_prime_guess, omega, ...
  r, R, V0, theta, B, cd_constant, c, beta_correction, fake_zero, i_max )

a = a_guess;
a_prime = a_prime_guess;

for i = 1: i_max
  tan_phi = (1 + a)*V0 / ((1 - a_prime)*omega*r); 
  phi = atan(tan_phi);
  alpha = theta - phi*180/(pi);
  
  cl = 0.1 * alpha + 0.4; % load cl
  cd = cd_constant; % load ct
  
  cn = cl*cos(phi) - cd*sin(phi); % compute cn
  ct = cl*sin(phi) + cd*cos(phi);

  sigma = sigma_function(c, B, r);
  F = tip_loss(B, R, r, phi); 
  a_star = (1 + a)*sigma*cn/(4*F*sin(phi)^2);
  a_new = beta_correction*a_star + (1 - beta_correction)*a; % update a 
  
  a_prime_star = (1 - a_prime)*sigma*ct/(4*F*sin(phi)*cos(phi)); 
  a_prime_new = beta_correction*a_prime_star + (1 - beta_correction)*a_prime;
  % update a_prime
  
  if (abs(a - a_new) < fake_zero) && (abs(a_prime - a_prime_new) < fake_zero)
    break
  end
  
  % update a and a_prime
  a = a_new;
  a_prime = a_prime_new;
end
end