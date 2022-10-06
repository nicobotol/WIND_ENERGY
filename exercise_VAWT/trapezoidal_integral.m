function [res] =  trapezoidal_integral(t_start, t_stop, y_vect, x_vect)

delta_x = x_vect(2) - x_vect(1);
s_start = round(t_start / delta_x); % sample where start the integration
s_stop = round(t_stop /delta_x);
s_size = s_stop-s_start-1;
res = 0;

for i=2:s_size
  res = res + 0.5*(y_vect(i) + y_vect(i-1))*delta_x;
end

end