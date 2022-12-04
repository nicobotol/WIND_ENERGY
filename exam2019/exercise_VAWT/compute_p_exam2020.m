function [Vrel, alpha] = compute_p_exam2020(V0, R, omega, theta, Wx, Wy)

x = -R*sin(theta);
y = R*cos(theta);

Vrel_x = omega*y + V0 - Wx;
Vrel_y = -omega*x + Wy;

Vnorm = (V0 - Wx)*sin(theta) - Wy*cos(theta);
Vtan = (V0 - Wx)*cos(theta) + Wy*sin(theta) + omega*R;
Vrel_sq = Vrel_x^2 + Vrel_y^2;
Vrel = sqrt(Vrel_sq);
%Vrel_sq = Vnorm^2 + Vtan^2;
alpha = atan(Vnorm / Vtan); % ange of attack (rad)
alphadeg = alpha*180/pi;

end