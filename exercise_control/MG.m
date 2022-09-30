function [Mg] = MG(omega, MG_coefficient, sincronous_velocity)
% omega (rad/s)
% sincronous_velocity (rad/s)

if omega < sincronous_velocity
  Mg = 0;
else
  Mg = MG_coefficient*(omega - sincronous_velocity);
end
end