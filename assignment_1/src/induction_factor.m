function [a, a_prime, ct, cn, phi] = induction_factor(a, a_prime, R, r, lambda, beta, Theta_p, cd, cl, sigma, B)
%   Function to calculate a and a_prime

  phi = atan(((1 - a)*R) / ((1 + a_prime) * lambda * r)); % flow angle

  Theta = beta + Theta_p; % Total pitch angle
  alpha = phi - Theta;  % angle of attack

  cn = cl*cos(phi) + cd*sin(phi); 
  ct = cl*sin(phi) - cd*cos(phi); 
  
  F = tip_loss(B, R, r, phi); % Prandtl's tip loss correction 

  %cT = (1 - a)^2*cn*sigma / (sin(phi)^2);

  % better guess for a and a_prime
  if a <= 1/3 
    %cT = 4*a*F*(1 - a); 
    a = (sigma * cn * (1 - a)) / (4 * F * sin(phi)^2);
  else
    %cT = 4*a*F*(1 - 1/4*(5 - 3*a)*a);
    cT = (1 - a)^2*cn*sigma / (sin(phi)^2);
    beta_correction = 0.1;
    a_star = cT / (4*F*(1 - 1/4*(5 - 3*a)*a));
    a = beta_correction*a_star + (1 - beta_correction)*a;
  end
  
  a_prime = sigma*ct*(1 + a_prime) / (4*F*sin(phi)*cos(phi));

end