%% Wind Turbine Technology and Aerodynamics - Assignment 2

clear 
close all
clc

%% Load data of the problem
parameters

%% Question 1
omega_M_rpm = 0.3:0.2:omega_M_max; % (rpm)

omega_M = [0:1:9 omega_M_max]; % rpm for the tables
%omega_M = [0:1:10]*pi/30;

omega_M = omega_M_rpm*pi/30; % (rad/s)
rpm = 30/pi*omega_M;

Pavailable = 0.5*rho*A*(omega_M*R/lambda_opt).^3;

Pmech = 0.5*rho*A*cp_opt*(omega_M*R/lambda_opt).^3; % mechanical power from the wind (W)
[Ig, Vg, fg, omega_E, delta] = rotor_equation(Rs, Ls, pp, flux, ...
  omega_M, Pmech);

% fig_mech_power = figure('Position', get(0, 'Screensize'));
% plot(rpm, Pmech, 'LineWidth', line_width)
% hold on
% plot(rpm, Pavailable,'LineWidth', line_width )
% hold off
% xlabel('Rotational speed [rpm]')
% ylabel('Mechanical power [W]')
% legend('Extracted power', 'Available power', 'Location','northwest')
% grid on
% title()
% set(gca, 'FontSize', font_size);
% saveas(fig_mech_power, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
%   '\WIND_ENERGY\assignment2\figures\fig_mech_power.png'],'png');

fig_current = figure('Position', get(0, 'Screensize'));
plot(rpm, Ig,'LineWidth', line_width)
xlabel('Rotational speed [rpm]')
ylabel('Generator current Ig [A]')
grid on
title('Generator current')
set(gca, 'FontSize', font_size);
saveas(fig_current, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\assignment2\figures\fig_current.png'],'png');

fig_voltage = figure('Position', get(0, 'Screensize'));
plot(rpm, Vg,'LineWidth', line_width)
xlabel('Rotational speed [rpm]')
ylabel('Generator voltage Vg [V]')
grid on
title('Generator voltage')
set(gca, 'FontSize', font_size);
saveas(fig_voltage, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\assignment2\figures\fig_voltage.png'],'png');

fig_frequency = figure('Position', get(0, 'Screensize'));
plot(rpm, fg,'LineWidth', line_width)
xlabel('Rotational speed [rpm]')
ylabel('Generator frequency [Hz]')
title('Generator frequency')
grid on
set(gca, 'FontSize', font_size);
saveas(fig_frequency, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\assignment2\figures\fig_frequency.png'],'png');


%% Question 2

% Losses
[P_loss, Q_loss, S_loss, P_g, Q_g, S_g] = losses(Ig, Vg, omega_E, Ls, Rs);

fig_mech_power = figure('Position', get(0, 'Screensize'));
plot(rpm, Pmech, 'LineWidth', line_width)
hold on
xlabel('Rotational speed [rpm]')
ylabel('Mechanical power [W]')
grid on
title('Mechanical power')
set(gca, 'FontSize', font_size);
saveas(fig_mech_power, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\assignment2\figures\fig_mech_power.png'],'png');


fig_losses_power = figure('Position', get(0, 'Screensize'));
plot(rpm, P_loss,'LineWidth', line_width)
hold on
plot(rpm, Q_loss,'LineWidth', line_width)
hold off
legend('Active [W]', 'Reactive [VAR]', 'Location','northwest')
xlabel('Rotational speed [rpm]')
ylabel('Power losses')
grid on
title('Power losses')
set(gca, 'FontSize', font_size);
saveas(fig_losses_power, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\assignment2\figures\fig_losses_power.png'],'png');


% Output complex power
fig_complex_power = figure('Position', get(0, 'Screensize'));
plot(rpm, S_g,'LineWidth', line_width)
hold on
plot(rpm, P_g,'LineWidth', line_width)
plot(rpm, Q_g,'LineWidth', line_width)
hold off
legend('Complex [VA]', 'Active [W]', 'Reactive [VAR]', 'Location','northwest')
xlabel('Rotational speed [rpm]')
ylabel('Power output')
title('Output power')
grid on
set(gca, 'FontSize', font_size);
saveas(fig_complex_power, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\assignment2\figures\fig_complex_power.png'],'png');

%% Question 3
clc
Vb = Vb_ll / sqrt(3); % phase voltage on VSC-B (V)

omega = 2*pi*f_grid; % angular velocity grid side

% [Ppoc, Qpoc] = transformer(Ib, Vb, Vpoc, Cc, Rc, Lc, omega, L2_prime, ...
%   R2_prime, n, Lm, R1, L1 );

Vpoc_prime_ll = Vpoc_ll*n ; % Vpoc
Vpoc_prime = Vpoc_prime_ll / sqrt(3); % transform the voltage to phase

Cc_prime = n^2*Cc;
Rc_prime = n^2*Rc;
Lc_prime = n^2*Lc;

% Z2_prime = (R2_prime + Rc_prime) + 1j*omega*(L2_prime + Lc_prime);
% Ipoc = 1/Z2_prime*(Vb - Ib*(R1 + 1j*omega*L1) - Vpoc_prime) - ...
%   1j*omega*Cc_prime*Vpoc_prime;
% Zm = 1j*omega*Lm;
% Z1 = R1+1j*omega*L1;
% Ipoc1 = Ib - (Vb-Ib*Z1)/Zm - 1j*Vpoc_prime*omega*Cc_prime;
Z1 = R1 + 1j*omega*L1;
Z_cable = Rc_prime + 1j*omega*Lc_prime;
Z2_prime = R2_prime +1j*L2_prime*omega;
Zm =  1j*omega*Lm;
Zc = 1/(1j*omega*Cc_prime);

for i =1:size(P_g,2)
  [Vb_sol(i), Vab_sol(i), Vcd_sol(i), Vpoc_sol(i), Ppoc(i), Qpoc(i)] = Q3(P_g(i), Z1, Z2_prime, Zm, Zc, Vpoc_prime, Z_cable);
end


% Pa = P_g(1) / 3; % phase power on the generator side [W]
% Pb = Pa;
% Qb = 0.2*Pb;
% % phi = atan(Qb/Pb); % leading angle
% % Sb_mag = sqrt(Pb.^2 + Qb.^2);
% Sb = Pb + 1j*Qb; % complex power on b side (VAR)
% threshold = 10;
% 
% for i = 1:1000
%   Vb = Vb_vector(i);
%   Ib = conj(Sb / Vb);
%   Vab = Vb - Ib*Z1;
%   Im = Vab/Zm;
%   I2 = Ib - Im;
%   Vpoc_star = Vab - I2*Z2_prime;
%   Vpoc_star_mag = abs(Vpoc_star);
%  
% 
%   if( abs(Vpoc_star_mag - Vpoc_prime)) < threshold
%     disp(strcat('convergence at', num2str(i)))
%     Vb_solution = Vb; 
%     Ic = Vpoc_star / Zc;
%     Ipoc = I2 - Ic;
%     Spoc = 3*Vpoc_star*conj(Ipoc)
%     Ppoc = real(Spoc)
%     Qpoc = imag(Spoc)
%     break
%   end
% 
% end


% Spoc = 3*Vpoc_prime.*conj(Ipoc); % the sign (-) comes from the fact
% % that current and voltage are assumed positive in opposite directions
% % Spoc = -3*Vpoc_prime_ll*conj(Ipoc);
% Ppoc = real(Spoc);
% Qpoc = imag(Spoc);

%%

% for p=1:47
% Pa = P_g(p);
% Pb = Pa/sqrt(3);
% Qb = 0.2*Pb;
% Sb = Pb + 1j*Qb;
% for i = 1:100
%   Vb = 4000/sqrt(3)*exp(1j*pi/180*i/100);
%   Ib = conj(Sb/Vb);
% 
%   Vab = Vb - Ib*Z1;
%   Im = Vab/Zm;
%   I2 = Ib - Im;
%   Vpoc_prime = Vab - I2*Z2_prime;
%   disp(abs(Vpoc_prime))
%   if abs(abs(Vpoc_prime) - Vpoc_ll/sqrt(3)*n) < 0.0000001 
%     disp('convergence ok')
%     break
%   else
%   end 
% 
% end
% Ic = Vpoc_prime ./ Zc;
% Ipoc = I2 - Ic;
% 
% Spoc = 3*Vpoc_prime.*conj(Ipoc); % the sign (-) comes from the fact
% Ppoc(p) = real(Spoc);
% Qpoc(p) = imag(Spoc);
% end

% POC active power
fig_POC_active = figure('Position', get(0, 'Screensize'));
plot(rpm, Ppoc,'LineWidth', line_width)
hold on
plot(rpm, Qpoc,'LineWidth', line_width)
hold off
legend('Active [W]', 'Reactive [VAR]', 'Location','northwest')
xlabel('Rotational speed [rpm]')
ylabel('Power at POC')
grid on
title('Power at POC')
set(gca, 'FontSize', font_size);
saveas(fig_POC_active, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\assignment2\figures\fig_POC_active.png'],'png');

% comparison of power
fig_POC_power = figure('Position', get(0, 'Screensize'));
plot(rpm, Ppoc,'LineWidth', line_width)
hold on
plot(rpm, P_g,'LineWidth', line_width)
plot(rpm, Pmech,'LineWidth', line_width)
hold off
legend('POC', 'Output', 'Mechanical', 'Location','northwest')
xlabel('Rotational speed [rpm]')
ylabel('Power [W]')
grid on
title('Power comparison')
set(gca, 'FontSize', font_size);
saveas(fig_POC_power, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\assignment2\figures\fig_POC_power.png'],'png');
%% Plot of voltages for question 3
% comparison of power
% fig_ = figure('Position', get(0, 'Screensize'));
figure()
plot(rpm, Vb_sol,'LineWidth', line_width)
hold on
plot(rpm, Vab_sol,'LineWidth', line_width)
plot(rpm, Vcd_sol,'LineWidth', line_width)
plot(rpm, Vpoc_sol,'LineWidth', line_width)
hold off
legend('vb', 'vab', 'vcd', 'Location','northwest')
xlabel('Rotational speed [rpm]')
ylabel('Voltage Vb')  
grid on


% set(gca, 'FontSize', font_size);
% saveas(fig_POC_power, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
%   '\WIND_ENERGY\assignment2\figures\fig_POC_power.png'],'png');
% 
% % Power for the farm of 8 WT
% Ppoc8 = 8*Ppoc;
% Qpoc8 = 8*Qpoc;
% fig_POC_active8 = figure('Position', get(0, 'Screensize'));
% plot(rpm, Ppoc8,'LineWidth', line_width)
% hold on
% plot(rpm, Qpoc8,'LineWidth', line_width)
% hold off
% legend('Active [W]', 'Reactive [VAR]', 'Location','northwest')
% xlabel('Rotational speed [rpm]')
% ylabel('Power at POC')
% grid on
% set(gca, 'FontSize', font_size);
% saveas(fig_POC_active8, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
%   '\WIND_ENERGY\assignment2\figures\fig_POC_active8.png'],'png');
% 
%% Question 4
% clc
% clear Ipoc
% syms Ipoc Ib Vb Ppoc I2 Vab Im Sb Pb
% % % phi = transformer2( Vb, Vpoc, Cc, Rc, Lc, omega, L2_prime, R2_prime, n, ...
% % %   Lm, R1, L1, Pb/3 );
% % % disp(num2str(phi))
% % 
% % % [Ppoc, Qpoc] = transformer(2.0165e+03-1.6132e+04i, 4000, Vpoc, Cc, Rc, Lc, omega, ...
% % %   L2_prime, R2_prime, n, Lm, R1, L1 );
% % % disp(Ppoc)
% % % disp(Qpoc)
% 
% Zm = 1j*omega*Lm;
% 
% Ppoc = Ipoc*Vpoc_prime;
% Ic = 1j*Vpoc_prime*omega*Cc_prime;
% I2 = Ic + Ipoc;
% Vab = Vpoc_prime + I2*(R2_prime + Rc_prime + 1j*(L2_prime + Lc_prime)*omega);
% Im = Vab/Zm;
% Ib = Im + I2;
% Vb = Vab + (R1 + 1j*omega*L1)*Ib;
% Sb = Vb*conj(Ib);
% Pb = real(Sb);
% 
% sol = vpasolve( Ppoc == Pb, Ipoc);
% 
% Ibsol = subs(Ib, Ipoc, sol);
% phi = eval(atan(imag(Ibsol)/real(Ibsol)))*180/pi
% 
% % figure()
% % plot(rpm, phi)
% 

clc
for i =1:size(P_g,2)
  [delta(i)] = Q4(P_g(i), Z1, Z2_prime, Zm, Zc, Vpoc_prime, Z_cable);
end

figure()
plot(rpm, delta)

%% Question 5
eta = Ppoc ./ Pmech;
%eta_tot = Ppoc'./ Pavailable;

fig_POC_power = figure('Position', get(0, 'Screensize'));
plot(rpm, eta,'LineWidth', line_width)
% hold on
% plot(rpm, eta_tot,'LineWidth', line_width)
% % hold off
% legend('\eta WT', '\eta total', 'Location', 'east')
xlabel('Rotational speed [rpm]')
ylabel('Efficiency \eta')
xlim([2 max(rpm)])
grid on
set(gca, 'FontSize', font_size);
saveas(fig_POC_power, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
  '\WIND_ENERGY\assignment2\figures\fig_POC_power.png'],'png');

%% Collect data in a table and plot it
% save_data(:, 1) = omega_M; % mechanical rotational speed [rad/s]
% save_data(:, 2) = omega_E; % electrical rotational speed [rad/s]
% save_data(:, 3) = Vg; % generator voltage [V]
% save_data(:, 4) = Ig; % generator current [A]
% save_data(:, 5) = fg; % genrator frequency [Hz]
% save_data(:, 6) = Pmech; % mechanical power [W]
% save_data(:, 7) = Pavailable; % available windspeed [m/s]
% save_data(:, 8) = P_loss; % active power loss in the generator [W]
% save_data(:, 9) = Q_loss; % reactive power loss in the generator [VAR]
% save_data(:, 10) = S_loss; % complex power loss in the generator [VA]
% save_data(:, 11) = S_g; % complex power output of the generator [VA]
% save_data(:, 12) = Ppoc; % Active power at the POC, for 1 WT [W]
% save_data(:, 13) = Qpoc; % Reactive power at the POC, for 1 WT [W]
% save_data(:, 14) = Ppoc8; % Active power at the POC, for 8 WT [W]
% save_data(:, 15) = Qpoc8; % Reactive power at the POC, for 8 WT [W]
% save_data(:, 16) = eta; % Efficiency of the WT (from mech power to POC)
% save_data(:, 17) = eta_tot; % Efficiency of the plant (from wind to POC)
% 
% columnLabels = {'\omega_M [rad/s]', '\omega_E [rad/s]', 'Vg [V]', 'Ig [A]'};
% 
% matrix2latex(save_data(:, 4), 'out.tex',...
%   'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', ...
%   'size', 'tiny');
