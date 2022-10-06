clear
close all
clc

%% Load airfoil data
load_data = importdata('airfoil.txt');
aoa = load_data.data(:, 1); % angle of attack
cl_data = load_data.data(:, 2); % cl vector
cd_data = load_data.data(:, 3); % cd vector

%% PARAMETERS
R = 3; % (m)
V0 = 8; % (m/s)
omega = 14; % (rad/s)
S = 0.2; % Solidity
rho = 1.225; % (kg/m^3)
B = 1;

blades = [1, 3]; % blades to investigate
blades_size = size(blades, 2); % number of blade configuration to investigate

delta_t = (5*pi)/(180*omega); % first guess fot time step increment
time_final = 10; % time when to stop the analysis
time_plot_initial = 0.0001; % time when start to plot values
time_plot_final = 9; % time when stop to plot the values

if (time_plot_final > time_final)
  error('Error in time length!')
end

ntime = ceil(time_final / delta_t); % number of time increment
delta_t = time_final / ntime; % time step increment

% time interval for plotting
s_initial = ceil(time_plot_initial/delta_t);
s_final = round(time_plot_final/delta_t);
s_plot = s_initial:s_final:1;

%%
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
  
  %     theta_i = 0; % clean theta_i for the next iteration
    end
  
    cT = px_sum / (rho*V0^2*R);
    cP = omega*pt_sum / (rho*V0^3);
    
    if a(n - 1) <= 1/3
      fg = 1;
    else
      fg = 0.25*(5 - a(n - 1));
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

% get the mean cP between 5 and 10 seconds
t_start = 5; % (s) time when start integration
t_stop = 10;
[sum_cp] =  trapezoidal_integral(t_start, t_stop, cp_vector, t);
cP_mean = sum_cp / (t_stop - t_start);

% plot the results
%build the legend
legend_name = strings(1, blades_size);
for b=1:blades_size
  legend_name(b) = strcat("B = ", num2str(blades(b)));
end

figure()
for b=1:blades_size
plot(t(s_initial:s_final), px_vector(s_initial:s_final, b));
hold on
end
hold off
legend(legend_name)
xlabel('Time (s)')
ylabel('px (N)')
title('px as function of time')

figure()
for b=1:blades_size
plot(t(s_initial:s_final), py_vector(s_initial:s_final, b));
hold on
end
hold off
legend(legend_name)
xlabel('Time (s)')
ylabel('py (N)')
title('py as function of time')

figure()
for b=1:blades_size
plot(t(s_initial:s_final), cp_vector(s_initial:s_final, b));
hold on
end
hold off
legend(legend_name)
xlabel('Time (s)')
ylabel('cP')
title('cP as function of time')

figure()
for b=1:blades_size
plot(t(s_initial:s_final), cT_vector(s_initial:s_final, b));
hold on
end
hold off
legend(legend_name)
xlabel('Time (s)')
ylabel('cT')
title('cT as function of time')

for b=1:blades_size
  figure()
  plot(t(s_initial:s_final), px_vector(s_initial:s_final, b));
  hold on
  plot(t(s_initial:s_final), py_vector(s_initial:s_final, b));
  hold off
  xlabel('Time (s)')
  ylabel('load (N/m)')
  legend('px', 'py')
  title(strcat('total load for B=', num2str(blades(b))))
end
