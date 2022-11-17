function [AEO, CF] = Annual_Energy_Output(P, V0, P_rated, A_weibull, k_weibull)
% compute the annual energy output

AEO = 0;  % Initialize AEO 

% compute the annual energy output
for i=1:(size(P, 2)-1)
  AEO = AEO + 0.5*(P(i) + P(i + 1)) * 8760 * ...
  ( exp(-(V0(i) / A_weibull)^k_weibull) - ...
  exp(-(V0(i + 1) / A_weibull)^k_weibull)  ); % (Wh)
end

% compute the capacity factor
CF = AEO / (8760 * P_rated);

end