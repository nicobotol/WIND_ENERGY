clear all 
close all
clc

% ALL THE VECTORS SHOULD BE ROW!

%% QUESTION 1

% % load the physical parameters from the file "parameters.m"
% parameters
% 
% % load the airfoil data from the different files
% [aoa_mat, cl_mat, cd_mat] = load_airfoil_data(filenames);
% 
% % load the blade data from "bladedat.txt"
% [r_vector, c_vector, beta_vetor, thick_vector] = load_blade_data(blade_filename);
% 
% % function to interpolate cl and cd for the given value of alpha and t/c
% [cl, cd] = interpolate_cl_cd(aoa_mat, cl_mat, cd_mat, thick_prof, deg2rad(-130), thick(18));
% 
% lambda = omega * R / V0;
% sigma = sigma_function(c, B, r); % sigma may be compute independently only once for each section r
% % compute the a and a_prime with the iterative method
% [a, a_prime, ct, cn, phi] = induction_factor_convergence(a, a_prime, R, r, lambda, beta, Theta_p, cd, cl, c, B, sigma, fake_zero, i_max);
% 
% V_rel = sqrt( (omega*r*(1 + a_prime))^2 + (V0*(1 - a))^2);
% 
% pn = 0.5 * rho * V_rel^2 * c * cn;
% pt = 0.5 * rho * V_rel^2 * c * ct;

% load the physical parameters from the file "parameters.m"
parameters

% load the airfoil data from the different files
[aoa_mat, cl_mat, cd_mat] = load_airfoil_data(filenames);

% load the blade data from "bladedat.txt"
[r_vector, c_vector, beta_vector, thick_vector] = load_blade_data(blade_filename);

r_item = size(r_vector, 2); % number of cross sections along the blade

lambda_vector = linspace(lambda_range(1), lambda_range(2), lambda_item); % vector of lambda equally distributed in the range
pitch_vector = linspace(pitch_range(1), pitch_range(2), pitch_item); % vector of pitch equally distributed in the range
cp = zeros(3, pitch_item*lambda_item); % matrix for store the couple (lambda, Theta) and the corresponding cp value

for l=1:lambda_item  % loop over lambda
  lambda = lambda_vector(l);
  for t=1:pitch_item % loop over pitch
    Theta_p = pitch_vector(t); % pitch angle
    clearvars a a_prime ct 
    cp_partial = zeros(1, r_item);
    %a = zeros(1, r_item);
    %a_prime = zeros(1, r_item);
    for i=1:r_item % loop over the blade positions  ATTENTION TO NOT BE TOO CLOSE TO TIP
      r = r_vector(i);
      beta = beta_vector(i);
      thick = thick_vector(i);
      c = c_vector(i);
      
      sigma = sigma_function(c, B, r);

      % compute the a and a_prime with the iterative method
      [a, a_prime, ct, ~, ~] = induction_factor_convergence(a_guess, a_prime_guess, R, r, lambda, beta, Theta_p, B, sigma, aoa_mat, cl_mat, cd_mat, thick_prof, thick, fake_zero, i_max);
      
      cp_partial(i) = r * ((1 - a)^2 + (lambda*r/R*(1 + a_prime))^2)*c*ct;
      
    end % stop looping over blade position r
    % integrate the partial results to get cp 
    pos = (l-1)*pitch_item + t;
%     a_prime
%     a
%     cp_partial
%     pause
    cp(3, pos) = lambda*B/A * trapezoidal_integral(r_vector, cp_partial);
    cp(1, pos) = lambda;
    cp(2, pos) = Theta_p;
  end % stop looping over pitch Theta_p
end % stop looping over lambda

%%
% plot the results
cp_vs_pitch = figure('Position', get(0, 'Screensize'));
legend_name_cp_vs_pitch = strings(1, lambda_item);
for l=1:lambda_item
  set = [1:1:pitch_item] + (l - 1)*pitch_item;
  plot(rad2deg(cp(2,set)), cp(3,set))
  legend_name_cp_vs_pitch(l) = strcat("\lambda = ", num2str(lambda_vector(l)));
  hold on
end
ylabel('cP')
xlabel('\Theta_P')
legend(legend_name_cp_vs_pitch)
title('cP as function of the pitch angle')
hold off

cp_vs_lambda = figure('Position', get(0, 'Screensize'));
legend_name_cp_vs_lambda = strings(1, pitch_item);
for p=1:pitch_item
  set2 = [1:pitch_item:size(cp,2)-pitch_item+1] + (p - 1);
  plot(cp(1,set2), cp(3,set2))
  legend_name_cp_vs_lambda(p) = strcat("\Theta = ", num2str(rad2deg(pitch_vector(p))));
  hold on
end
ylabel('cP')
xlabel('\lambda')
legend(legend_name_cp_vs_lambda)
title('cP as function of \lambda')
hold off



% procedue to create a matrix of cP as function of lamda and theta
%
% create a vector of lambda within the suggested range
% create a vector of theta within the suggested range
% loop over lambda
%   loop over theta
%     loop over r (not arrive to close to the tip)
%       compute sigma
%         loop until convergence to get a, a_prime, ct
%         end loop until convergence
%       store a, a_prime, ct
%     end loop over r
%     integrate and find the global cP
%     store the result in the column of a matrix
%   end loop over theta
%   store the results in the row of a matrix
% end loop over lambda

%% QUESTION 2

% optimal lambda
[cp_max, cp_max_pos] = max(cp(3,:));
lambda_opt = cp(1, cp_max_pos);
Theta_p_opt = cp(2, cp_max_pos);
% rated velocity
V0_rated = (P_rated / (0.5*cp_max/100*rho*A) )^(1/3); % rated wind velocity (m/s)
omega_max = V0_rated * lambda_opt / R;

V0_step = (V0_rated - V0_cutin) / (V0_item - 1);
V0_vector = [V0_cutin: V0_step: V0_rated];
%V0_vector = linespace(V0_cutin, V0_rated, V0_item);
omega_plot = lambda_opt / R * V0_vector;

omega_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_vector, omega_plot);
xlabel("Wind velocity (m/s)")
ylabel("\omega (rad/s)")
title("Rotational speed as function of wind velocity")

%% QUESTION 3
% first part of the question
Theta_p = Theta_p_opt;

% second part of the question 
P = size(1, V0_item); % initialize vector of power
T = size(1, V0_item); % initialize vector of thrust
cP = size(1, V0_item);
cT = size(1, V0_item);

V0_vector_cut_in_out = linspace(V0_cutin, V0_cut_out, V0_cut_in_out_item);
for v=1:V0_cut_in_out_item % loop over differnet velocities
  V0_actual = V0_vector_cut_in_out(v);
  
  pn_partial = zeros(1, r_item); % initialize the pn vector
  pt_partial = zeros(1, r_item); % initialize the pn vector

  for i=1:r_item % loop over the blade positions  ATTENTION TO NOT BE TOO CLOSE TO TIP
    r = r_vector(i);
    beta = beta_vector(i);
    thick = thick_vector(i);
    c = c_vector(i);
    sigma = sigma_function(c, B, r);
    
    if V0_actual < V0_rated
      lambda = lambda_opt;
      omega_actual = V0_actual * lambda / R;
    else
      omega_actual = omega_max;
      lambda = omega_actual * R / V0_actual;
    end
    %
    % I have not choose which Theta_p to use, it must be done accordingly
    % to the first part of this question
    %
    % compute the a and a_prime with the iterative method
    [a, a_prime, ct, cn, ~] = induction_factor_convergence(a_guess, a_prime_guess, R, r, lambda, beta, Theta_p, B, sigma, aoa_mat, cl_mat, cd_mat, thick_prof, thick, fake_zero, i_max);
    
    Vrel = sqrt((V0_actual*(1 - a))^2 + (omega_actual*r*(1 + a_prime))^2); % (m/s)

    pn_partial(i) = 0.5*rho*c*cn*Vrel^2;
    pt_partial(i) = r * 0.5*rho*c*ct*Vrel^2;
  end % stop looping over blade position r

  T(v) = B*trapezoidal_integral(r_vector, pn_partial); % thrust (N)
  P(v) = omega_actual*B*trapezoidal_integral(r_vector, pt_partial); % power (W)

  cT(v) = T(v) / (0.5*rho*V0_actual^2*A);
  cP(v) = P(v) / (0.5*rho*V0_actual^3*A);
end % stop looping over different velocities

P_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_vector_cut_in_out, P);
xlabel('Wind velocity V0 (m/s)')
ylabel('Power (W)')

T_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_vector_cut_in_out, T);
xlabel('Wind velocity V0 (m/s)')
ylabel('Thrust (N)')

cP_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_vector_cut_in_out, cP);
xlabel('Wind velocity V0 (m/s)')
ylabel('cP')

cT_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_vector_cut_in_out, cT);
xlabel('Wind velocity V0 (m/s)')
ylabel('cT')