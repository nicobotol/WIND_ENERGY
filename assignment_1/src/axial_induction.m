function [a, a_prime] = induction_factor(a, a_prime, R, r, lambda, beta, Theta_p, cd, cl, sigma)
%   Function to calculate a and a_prime

  phi = atan((1 - a) / (1 + a_prime) * R / r / lambda); % flow angle

  Theta = beta + Theta_p; % Total pitch angle
  alpha = phi - Theta;  % angle of attack

  cn = cl*cos(phi) + cd*sin(phi); 
  ct = cl*sin(phi) - cd*cos(phi); 
  
  F = tip_loss(B, R, r, phi); % Prandtl's tip loss correction 

  cT = (1 - a)^2*cn*sigma / sin(phi)^2;

  % better guess for a and a_prime
  if a <= 0.33 
    a = (sigma * cn * (1 - a)) / (4 * F * sin(phi)^2);
  else
    beta_correction = 0.1;
    a_star = cT / (4*F*(1 - 0.25*(5 - 3*a)*a));
    beta_correction*a_star + (1 - beta_correction)*a;
  end
  
  a_prime = sigma*ct*(1 + a_prime) / (4*F*sin(phi)*cos(phi));

end