function [F] = tip_loss(B, R, r, phi)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
  F = 2/pi * acos( exp(- B * (R - r) / (2*r*sin(abs(phi)))));
end