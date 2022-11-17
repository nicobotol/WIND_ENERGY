function [h] = weibull_pdf(V0, A_weibull, k_weibull)
% probability density function of weibull distribution

h = k_weibull/A_weibull*(V0/A_weibull)^(k_weibull - 1)*exp(-(V0/A_weibull)^k_weibull);
end