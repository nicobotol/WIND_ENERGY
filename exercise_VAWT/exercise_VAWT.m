clear
close all
clc

%% Load airfoil data
load_data = importdata('airfoil.txt');
aoa = load_data.data(:, 1); % angle of attack
cl_data = load_data.data(:, 2); % cl vector
cd_data = load_data.data(:, 3); % cd vector

font_size = 15; 
line_size = 1.75;

%% HELPER LINES FOR PRINTING
% fig_Mg = figure('Position', get(0, 'Screensize'));
% set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
% saveas(fig_Mg, ['C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO' ...
%   '\WIND_ENERGY\exercise_control\images\fig_Mg.png'],'png');


%% PARAMETERS
R = 3; % (m)
V0 = 8; % (m/s)
omega = 14; % (rad/s)
S = 0.2; % Solidity
rho = 1.225; % (kg/m^3)

blades = [1 2 3]; % blades to investigate
blades_size = size(blades, 2); % number of blade configuration to investigate
delta_t = (1*pi)/(180*omega); % first guess fot time step increment
time_final = 20; % time when to stop the analysis
time_plot_initial = 0.0006; % time when start to plot values
time_plot_final = 9; % time when stop to plot the values
t_start = 4.5; % (s) time when start integration
t_stop = 20;  % (s) time when stop integration
if (time_plot_final > time_final)
  error('Error in time length!')
end

ntime = ceil(time_final / delta_t); % number of time increment
delta_t = time_final / ntime; % time step increment

% time interval for plotting
s_initial = ceil(time_plot_initial/delta_t);
s_final = round(time_plot_final/delta_t);
s_plot = s_initial:s_final:1;

%% Do all the computations
a = zeros(ntime, 1);
theta1 = zeros(ntime, 1);
t = zeros(ntime, 1);
px_vector = zeros(ntime, blades_size);
py_vector = zeros(ntime, blades_size);
cp_vector = zeros(ntime, blades_size);
cT_vector = zeros(ntime, blades_size);
px_sum = 0;
py_sum = 0;
pt_sum = 0;
% struct_loads{} = zeros(ntime, 2*blades_size);

% initialize a, theta1, t
a(1) = 0;
theta1(1) = 0;
t(1) = 0;

% compute the time constant
tau = 2*R/V0;

for b=1:blades_size
  B = blades(b); %select the numbe rof blades
  
  for n=2:ntime % run over time
    t(n) = (n - 1)*delta_t;
    theta1(n) = theta1(n - 1) + omega*delta_t;
    W_x0 = a(n - 1)*V0;
    
    % clear variables from previous cycle
    px_sum = 0;
    py_sum = 0;
    pt_sum = 0;
  
    for i=1:B % run over number of blades
      c = S*R/B;
      theta_i = theta1(n) + 2*pi*(i - 1)/B;
      [Wx, Wy] = compute_W(theta_i, W_x0);
      [px, py, pt] = compute_p(V0, R, rho, omega, theta_i, c, aoa, cl_data, cd_data, Wx, Wy);
      
      px_sum = px_sum + px;
      py_sum = py_sum + py;
      pt_sum = pt_sum + pt;

      % save load blade by blade
      struct_loads{b}.p(n, 3*i - 2) = px; 
      struct_loads{b}.p(n, 3*i - 1 ) = py;
      struct_loads{b}.p(n, 3*i ) = pt;

  %     theta_i = 0; % clean theta_i for the next iteration
    end
  
    cT = px_sum / (rho*V0^2*R);
    cP = omega*pt_sum / (rho*V0^3);
    
    if a(n - 1) <= 1/3
      fg = 1;
    else
      fg = 0.25*(5 - 3*a(n - 1));
    end
  
    aqs = cT / (4*(1 - fg*a(n - 1)));
    a(n) = aqs + (a(n - 1) - aqs)*exp(-delta_t / tau);
  
    % save the data
    px_vector(n, b) = px_sum;
    py_vector(n, b) = py_sum;
    cp_vector(n, b) = cP;
    cT_vector(n, b) = cT;
  
  end
end


%% get the mean cP between 5 and 10 seconds

for i = 1:blades_size
  [sum_cp] =  trapezoidal_integral(t_start, t_stop, cp_vector(:, i), t);
  [sum_cT] =  trapezoidal_integral(t_start, t_stop, cT_vector(:, i), t);
  cP_mean = sum_cp / (t_stop - t_start);
  cT_mean = sum_cT / (t_stop - t_start);
  disp( strcat('For B=', num2str(blades(i)), ': cp=', num2str(cP_mean), '; cT=', num2str(cT_mean)))
end
cd
for i = 1:blades_size
  [sum_cp] =  trapezoidal_integral(t_stop - 2*pi/omega, t_stop, cp_vector(:, i), t);
  [sum_cT] =  trapezoidal_integral(t_stop - 2*pi/omega, t_stop, cT_vector(:, i), t);
  cP_mean = sum_cp / (2*pi/omega);
  cT_mean = sum_cT / (2*pi/omega);
  disp( strcat('For B=', num2str(blades(i)), ': cp=', num2str(cP_mean), '; cT=', num2str(cT_mean)))
end
%% plot the results

% build the legend
legend_name = strings(1, blades_size);
for b=1:blades_size
  legend_name(b) = strcat("B = ", num2str(blades(b)));
end

%%
% plot px for every blade number
figure()
for b=1:blades_size
plot(t(s_initial:s_final), px_vector(s_initial:s_final, b), 'LineWidth', line_size);
hold on
end
hold off
legend(legend_name, 'FontSize', font_size)
xlabel('Time (s)')
ylabel('px (N)')
title('px as function of time')
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
grid on

%%
% plot py for every blade number
figure()
for b=1:blades_size
plot(t(s_initial:s_final), py_vector(s_initial:s_final, b), 'LineWidth', line_size);
hold on
end
hold off
legend(legend_name, 'FontSize', font_size)
xlabel('Time (s)')
ylabel('py (N)')
title('py as function of time')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)

%%
% plot px and py of one blade on the same graph
for b=1:blades_size
  figure()
  plot(t(s_initial:s_final), px_vector(s_initial:s_final, b), 'LineWidth', line_size);
  hold on
  plot(t(s_initial:s_final), py_vector(s_initial:s_final, b), 'LineWidth', line_size);
  hold off
  xlabel('Time (s)')
  ylabel('load (N/m)')
  legend('px', 'py', 'FontSize', font_size)
  title(strcat('total load for B=', num2str(blades(b))))
  grid on
  set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
end

%%
% plot cp for every blade
figure()
for b=1:blades_size
plot(t(s_initial:s_final), cp_vector(s_initial:s_final, b), 'LineWidth', line_size);
hold on
end
hold off
legend(legend_name, 'FontSize', font_size)
xlabel('Time (s)')
ylabel('cP (-)')
title('cP as function of time')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)

%%
% plot cT for every blade
figure()
for b=1:blades_size
plot(t(s_initial:s_final), cT_vector(s_initial:s_final, b), 'LineWidth', line_size);
hold on
end
hold off
legend(legend_name, 'FontSize', font_size)
xlabel('Time (s)')
ylabel('cT (-)')
title('cT as function of time')
grid on
set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)

%%
% plot cp and cT of one blade on the same graph
for b=1:blades_size
  figure()
  plot(t(s_initial:s_final), cp_vector(s_initial:s_final, b), 'LineWidth', line_size);
  hold on
  plot(t(s_initial:s_final), cT_vector(s_initial:s_final, b), 'LineWidth', line_size);
  hold off
  xlabel('Time (s)')
  ylabel('coeff. (-)')
  legend('c_P', 'c_T', 'FontSize', font_size)
  title(strcat('c_P and c_T for B=', num2str(blades(b))))
  grid on
  set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
end

%%
% plot px on all the blades for one blade configuration
for b=1:blades_size
  figure()
  for i=1:blades(b)
    plot(t(s_initial:s_final), struct_loads{b}.p((s_initial:s_final), i*3 - 2), 'LineWidth', line_size);
    hold on
  end
  hold off
  xlabel('Time (s)')
  ylabel('load (N/m)')
  legend(legend_name, 'FontSize', font_size)
  title(strcat('px load for B=', num2str(blades(b))))
  grid on
  set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
end

%%
% plot py on all the blades for one blade configuration
for b=1:blades_size
  figure()
  for i=1:blades(b)
    plot(t(s_initial:s_final), struct_loads{b}.p((s_initial:s_final), i*3 - 1), 'LineWidth', line_size);
    hold on
  end
  hold off
  xlabel('Time (s)')
  ylabel('load (N/m)')
  legend(legend_name, 'FontSize', font_size)
  title(strcat('py load for B=', num2str(blades(b))))
  grid on
  set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
end

%%
% plot pt on all the blades for one blade configuration
for b=1:blades_size
  figure()
  for i=1:blades(b)
    plot(t(s_initial:s_final), struct_loads{b}.p((s_initial:s_final), i*3), 'LineWidth', line_size);
    hold on
  end
  hold off
  xlabel('Time (s)')
  ylabel('load (N/m)')
  legend(legend_name, 'FontSize', font_size)
  title(strcat('pt load for B=', num2str(blades(b))))
  grid on
  set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
end

%%
% plot pt on all the blades for one blade configuration
for b=1:blades_size
  figure()
  for i=1:blades(b)
    plot(t(s_initial:s_final), struct_loads{b}.p((s_initial:s_final), i*3), 'LineWidth', line_size);
    hold on
  end
  hold off
  xlabel('Time (s)')
  ylabel('load (N/m)')
  legend(legend_name, 'FontSize', font_size)
  title(strcat('pt load for B=', num2str(blades(b))))
  grid on
  set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
end

%%
% plot pt on 1 blade for one blade configuration
% plot pt on all the blades for one blade configuration
for b=1:blades_size
  figure()
  plot(t(s_initial:s_final), struct_loads{b}.p((s_initial:s_final), 3), 'LineWidth', line_size);
  hold on
  hold off
  xlabel('Time (s)')
  ylabel('load (N/m)')
  legend(legend_name, 'FontSize', font_size)
  title(strcat('pt load on one blade, for B=', num2str(blades(b))))
  grid on
  set(gca, 'FontAngle', 'oblique', 'FontSize', font_size)
end
