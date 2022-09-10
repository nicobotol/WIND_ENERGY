function [out] = trapezoidal_integral(x_vect, y_vect)
% This function implements the integration with the trapezoidal rules

out = 0;
for i = 1:size(x_vect)
  out = out + ((y_vect(i) + y_vect(i + 1))*(x_vect(i + 1) - x_vect(i))) / 2;
end