function [cP] = cP_section(B, c, lambda, ct, a, R, phi)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

cP = B*c*lambda*ct*(1 - a)^2 / (2*pi*R*sin(phi)^2);
end