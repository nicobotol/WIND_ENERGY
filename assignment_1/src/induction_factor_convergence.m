function [a, a_prime, ct, cn, phi] = induction_factor_convergence(a, a_prime, R, r, lambda, beta, Theta_p, cd, cl, c, B, fake_zero, i_max)
% This function performs the calculation to make the induction functions
% converge

for i = 1:i_max
  % initialize a and a_prime
  a_old = a;
  a_prime_old = a_prime;
  
  % update a and a_prime
  sigma = sigma_function(c, B, r);
  [a, a_prime, ct, cn, phi] = induction_factor(a, a_prime, R, r, lambda, beta, Theta_p, cd, cl, sigma, B);
  
  % compute the error
  epsilon = abs(a_old - a);
  epsilon_prime = abs(a_prime_old - a_prime);
  
  % brake the code if convergence
  if (epsilon < fake_zero) && (epsilon_prime < fake_zero)
    break
  end

end

end