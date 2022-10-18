function [sum_T, sum_M] = trapezoidal_integral(dT, dM, r_vector, r_size)
  sum_T = 0;
  sum_M = 0;

  for i = 2:r_size
    sum_T = sum_T + (dT(i - 1) + dT(i)) * (r_vector(i) - r_vector(i - 1))*0.5;
    sum_M = sum_M + (dM(i - 1) + dM(i)) * (r_vector(i) - r_vector(i - 1))*0.5;
  end

end