% This code computes the Annual Energy Production provided the power
% produced in a range of velocities, the velocities and the parameters of
% the Weibull distribution.
% Furthermore, if the rated power is provided, then the capacity facotr
% could be computed

%%
clear 
close all
clc

%% Defintion of parameters

% Linewidth and tett size for plot
line_width = 2;
font_size = 10;

% parameters of the weibull distribution
A_weibull = 8;
k_weibull = 2;

%rated power
P_rated = 10.64e6; % (W)

% wind speed and power
V0_start = 5; % velocity where to start the weibull
V0_stop = 7.5; % velocity where to stop the weibull 
P0 = 2e6; % power in the range of velocity (W)

%%
V0 = V0_start:0.1:V0_stop; % build the vector with wind velocity
P = P0*ones(1, size(V0, 2)); % build the vector with the power of each 
% wind velocity (W)
[AEO, CF] = Annual_Energy_Output(P, V0, P_rated, A_weibull, k_weibull);

% In case of multiple power for multiple wind speed
% V1 = 10:0.1:13; % build the vector with wind velocity
% P1 = 2e6*ones(1, size(V1, 2)); % build the vector with the power of each 
% % wind velocity (W)
% [AEO1, CF1] = Annual_Energy_Output(P1, V1, P_rated, A_weibull, k_weibull);
% 
% V2 = 13:0.1:25; % build the vector with wind velocity
% P2 = 3e6*ones(1, size(V2, 2)); % build the vector with the power of each 
% % wind velocity (W)
% [AEO2, CF2] = Annual_Energy_Output(P2, V2, P_rated, A_weibull, k_weibull);
% AEO = AEO + AEO1 + AEO2; % sum the energy production in the three part

disp(strcat('The total annula energy production is ', num2str(AEO/1e9), ...
  ' GWh'));

% build and plot the weibull distribution
weibull = zeros(1, size(V0, 2)); % initialize the vector to store the pdf
for i=1:size(V0, 2)
  weibull(i) = weibull_pdf(V0(i), A_weibull, k_weibull);
end

% print the weibull distribution
weibull_plot = figure('Position', get(0, 'Screensize'));
plot(V0, weibull, 'LineWidth', line_width) 
xlabel('Wind velocity V_0 (m/s)')
ylabel('Weibull PDF')
title('Weibull pdf as function of the wind velocity')
ax = gca;
ax.FontSize = font_size;
saveas(weibull_plot, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\' ...
  'WIND_ENERGY\exercise_annual_energy_production\weibull_plot.png'],'png');