function [Mr] = MR(rho, V0, R, cP, omega)

Mr = 0.5*rho*V0^3*pi*R^2*cP / omega;

end