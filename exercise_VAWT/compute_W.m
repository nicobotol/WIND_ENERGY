function [Wx, Wy] = compute_W(theta_i, W_x0)
  Wx = W_x0*(1 - 0.4*sin(theta_i));
  Wy = 0.4*Wx*cos(theta_i);
end