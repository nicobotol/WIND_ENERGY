function [px, py, pt] = compute_p(V0, R, rho, omega, theta, c, aoa, cl_data, cd_data, Wx, Wy)

x = -R*sin(theta);
y = R*cos(theta);

Vrel_x = omega*y + V0 - Wx;
Vrel_y = -omega*x + Wy;

Vnorm = (V0 - Wx)*sin(theta) - Wy*cos(theta);
Vtan = (V0 - Wx)*cos(theta) + Wy*sin(theta) + omega*R;
Vrel_sq = Vrel_x^2 + Vrel_y^2;
%Vrel_sq = Vnorm^2 + Vtan^2;
alpha = atan(Vnorm / Vtan); % (rad)
alphadeg = alpha*180/pi;

% cl = interp1(deg2rad(aoa), cl_data, alpha);
% cd = interp1(deg2rad(aoa), cd_data, alpha);

cl = interp1(aoa, cl_data, alphadeg);
cd = interp1(aoa, cd_data, alphadeg);

l = 0.5*rho*Vrel_sq*c*cl;
d = 0.5*rho*Vrel_sq*c*cd;

cos_beta = Vrel_y / sqrt(Vrel_sq);
sin_beta = Vrel_x / sqrt(Vrel_sq); 

px = l*cos_beta + d*sin_beta;
py = -l*sin_beta + d*cos_beta;

% phi = alpha + theta;
% pt = l*sin(alpha) - d*cos(alpha);

pt = -(px*cos(theta) + py*sin(theta)); 
end