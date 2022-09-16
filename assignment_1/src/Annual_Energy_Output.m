function [AEO_25, AEO_20, CF_25, CF_20, delta_AOE] = Annual_Energy_Output(P, V0, P_rated, A_weibull, k_weibull)
% compute the annual energy output

AEO_25 = 0;  % Initialize AEO
AEO_20 = 0; % 

for i=1:(size(P, 2)-1)
  AEO_25 = AEO_25 + 0.5*(P(i) + P(i + 1)) * 8760 * (  exp(-(V0(i) / A_weibull)^k_weibull) - exp(-(V0(i + 1) / A_weibull)^k_weibull)  ); % (Wh)
end

CF_25 = AEO_25 / (8760 * P_rated);

% second part: energy output up to V0 = 20 (m/s)
index_20 = find(V0 == 20);
P_20 = P(1:index_20);
V0_20 = V0(1:index_20);

for i=1:(size(P_20, 2)-1)
  AEO_20 = AEO_20 + 0.5*(P_20(i) + P_20(i + 1)) * 8760 * (  exp(-(V0_20(i) / A_weibull)^k_weibull) - exp(-(V0_20(i + 1) / A_weibull)^k_weibull) ); % (W/h)
end

CF_20 = AEO_20 / (8760 * P_rated);

delta_AOE = AEO_25 - AEO_20;

end