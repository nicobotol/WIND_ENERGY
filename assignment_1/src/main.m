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
r_item_no_tip = r_item - 1; % numer of cross section to take into account

lambda_vector = linspace(lambda_range(1), lambda_range(2), lambda_item); % vector of lambda equally distributed in the range
pitch_vector = linspace(pitch_range(1), pitch_range(2), pitch_item); % vector of pitch equally distributed in the range
cP_cT_mat = zeros(4, pitch_item*lambda_item); % matrix for store the couple (lambda, Theta) and the corresponding cp and cn values

contour_mat_cP = zeros(pitch_item, lambda_item); % matrix of cPs for the contour plot
contour_mat_cT = zeros(pitch_item, lambda_item); % matrix of cTs for the contour plot

for l=1:lambda_item  % loop over lambda
  lambda = lambda_vector(l);
  for t=1:pitch_item % loop over pitch
    Theta_p = pitch_vector(t); % local pitch angle
    clearvars a a_prime ct 
    cp_partial = zeros(1, r_item_no_tip);
    cT_partial = zeros(1, r_item_no_tip);

    [cp_partial, cT_partial,~,~] = cP_cT_partial(r_item_no_tip, r_vector, ...
    beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, R, lambda, ...
    Theta_p, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero,rho, 0,0, i_max);

    % integrate the partial results to get cp 
    pos = (l-1)*pitch_item + t;
    
    cP = lambda*B/(A*R) * trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial); 
    cT = B / A * trapezoidal_integral(r_vector(1:r_item_no_tip), cT_partial);
    cP_cT_mat(3, pos) = cP; %cP
    cP_cT_mat(4, pos) = cT; %cT
    cP_cT_mat(1, pos) = lambda;
    cP_cT_mat(2, pos) = Theta_p;

    contour_mat_cP(t, l) = cP;
    contour_mat_cT(t,l) = cT;
  end % stop looping over pitch Theta_p
end % stop looping over lambda

% optimal lambda
[cp_max, cp_max_pos] = max(cP_cT_mat(3,:));
lambda_opt = cP_cT_mat(1, cp_max_pos);
Theta_p_opt = cP_cT_mat(2, cp_max_pos);

disp_lambda = strcat('\lambda corresponding to max cP is \lambda=', num2str(lambda_opt));
disp_Theta_p = strcat('Pitch corresponding to max cP is \Theta_p=', num2str(rad2deg(Theta_p_opt)));
disp(disp_lambda)
disp(disp_Theta_p)
%%

% plot the results
legend_name_cp_vs_pitch = strings(1, lambda_item);
cp_vs_pitch = figure('Position', get(0, 'Screensize'));
for l=1:lambda_item
  set1 = [1:1:pitch_item] + (l - 1)*pitch_item;
  plot(rad2deg(cP_cT_mat(2,set1)), cP_cT_mat(3,set1))
  legend_name_cp_vs_pitch(l) = strcat("\lambda = ", num2str(lambda_vector(l)));
  hold on
end
ylabel('cP')
xlabel('\Theta_P')
legend(legend_name_cp_vs_pitch)
title('cP as function of the pitch angle')
hold off
ax = gca;
ax.FontSize = font_size;
saveas(cp_vs_pitch, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\cp_vs_pitch.png','png');
%%
cp_vs_lambda = figure('Position', get(0, 'Screensize'));
legend_name_cp_vs_lambda = strings(1, pitch_item);
for p=1:pitch_item
  set2 = [1:pitch_item:size(cP_cT_mat,2)-pitch_item+1] + (p - 1);
  plot(cP_cT_mat(1,set2), cP_cT_mat(3,set2))
  legend_name_cp_vs_lambda(p) = strcat("\Theta = ", num2str(rad2deg(pitch_vector(p))),"°");
  hold on
end
ylabel('cP')
xlabel('\lambda')
legend(legend_name_cp_vs_lambda)
title('cP as function of \lambda')
hold off
ax = gca;
ax.FontSize = font_size;
saveas(cp_vs_lambda, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\cp_vs_lambda.png','png');

% plot the results
cT_vs_pitch = figure('Position', get(0, 'Screensize'));
legend_name_cp_vs_pitch = strings(1, lambda_item);
for l=1:lambda_item
  set3 = [1:1:pitch_item] + (l - 1)*pitch_item;
  plot(rad2deg(cP_cT_mat(2,set3)), cP_cT_mat(4,set3))
  legend_name_cp_vs_pitch(l) = strcat("\lambda = ", num2str(lambda_vector(l)));
  hold on
end
ylabel('cT')
xlabel('\Theta_P')
legend(legend_name_cp_vs_pitch)
title('cT as function of the pitch angle')
hold off
ax = gca;
ax.FontSize = font_size;
saveas(cT_vs_pitch, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\cT_vs_pitch.png','png');

cT_vs_lambda = figure('Position', get(0, 'Screensize'));
legend_name_cp_vs_lambda = strings(1, pitch_item);
for p=1:pitch_item
  set4 = [1:pitch_item:size(cP_cT_mat,2)-pitch_item+1] + (p - 1);
  plot(cP_cT_mat(1,set4), cP_cT_mat(4,set4))
  legend_name_cp_vs_lambda(p) = strcat("\Theta = ", num2str(rad2deg(pitch_vector(p))), "°");
  hold on
end
ylabel('cT')
xlabel('\lambda')
legend(legend_name_cp_vs_lambda)
title('cT as function of \lambda')
hold off
ax = gca;
ax.FontSize = font_size;
saveas(cT_vs_lambda, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\cT_vs_lambda.png','png');

% contour plot
contour_plot_cP = figure('Position', get(0, 'Screensize'));
contourf(lambda_vector, rad2deg(pitch_vector), contour_mat_cP); % To display value on the plot use ,'ShowText','on'
hold on
plot(lambda_opt, Theta_p_opt, 'r.', 'MarkerSize',30)
text(7.9, 0, num2str(cp_max), 'Color','r', 'FontSize', font_size)
hold off
colorbar()
xlabel('\lambda')
ylabel('\Theta_p (°)')
title('Contour plot of c_P')
ax = gca;
ax.FontSize = font_size;
saveas(contour_plot_cP, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\contour_plot_cP.png','png');

contour_plot_cT = figure('Position', get(0, 'Screensize'));
contourf(lambda_vector, rad2deg(pitch_vector), contour_mat_cT);
colorbar()
xlabel('\lambda')
ylabel('\Theta_p (°)')
title('Contour plot of c_T')
ax = gca;
ax.FontSize = font_size;
saveas(contour_plot_cT, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\contour_plot_cT.png','png');

disp('Ex 1 done')
%% QUESTION 2

% rated velocity
V0_rated = (P_rated / (0.5*cp_max*rho*A) )^(1/3); % rated wind velocity (m/s)
omega_max = V0_rated * lambda_opt / R;

disp_V = strcat('The rated velocity V0 is V0=', num2str(V0_rated), ' (m/s)' );
disp_omega_max = strcat('The maximum rotational speed is \omega=', num2str(omega_max), ' (rad/s)' );
disp(disp_V)
disp(disp_omega_max)

V0_step = (V0_rated - V0_cutin) / (V0_item - 1);
V0_vector = [V0_cutin: V0_step: V0_rated];
%V0_vector = linespace(V0_cutin, V0_rated, V0_item);
omega_plot = lambda_opt / R * V0_vector;
%%
omega_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_vector, omega_plot, 'LineWidth', line_width);
hold on
plot([V0_rated V0_cut_out], [omega_max omega_max], 'LineWidth', line_width)
hold off
legend('V0 < V0_{rated}', 'V0 > V0_{rated}', 'Location','southeast' )
xlabel("Wind velocity (m/s)")
ylabel("\omega (rad/s)")
title("Rotational speed as function of wind velocity")
ylim([0.3 1.1])
ax = gca;
ax.FontSize = font_size;
saveas(omega_vs_V0, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\omega_vs_V0.png','png');



% %% QUESTION 3
% % first part of the question
% 
% % Method 1: Interpolate the 3d plot and look for the data
% % F = scatteredInterpolant(cP_cT_mat(3,:)', cP_cT_mat(1,:)', cP_cT_mat(2,:)'); %(cP, lambda, Theta_p)
% % rad2deg(F(0.42, 9.4))
% 
% % Method 2: look in the table of values
% % V0 = 20;
% % lambda = omega_max * R / V0;
% % cP = P_rated / (0.5*rho*V0^3*A);
% 
% % Method 3: interpolate with a polynomial
% V0_vector_cut_in_out = linspace(V0_cutin, V0_cut_out, V0_cut_in_out_item);
%                                                                                                          
% pitch_vector_e3 = linspace(pitch_range_e3(1), pitch_range_e3(2), pitch_item_e3); % vector of pitch equally distributed in the range 
% P_e3 = zeros(1, pitch_item_e3);
% 
%   V0_e3 = 20;
%   % compute lambda
%   if V0_e3 < V0_rated
%     lambda = lambda_opt;
%     omega_e3 = V0_e3 * lambda / R;
%   else
%     omega_e3 = omega_max;
%     lambda = omega_e3 * R / V0_e3;
%   end
% 
%   % chose Theta_p
%   
%   for p = 1:pitch_item_e3 %loop over different pitch
%     Theta_p_e3 = pitch_vector_e3(p);
% 
%     [cp_partial, cT_partial] = cP_cT_partial(r_item_no_tip, r_vector, ...
%     beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, R, lambda, ...
%     Theta_p_e3, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, i_max);
% 
%     P_e3(p) = 0.5*omega_e3*B*rho*V0_e3^2*trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial);
%     
%   end
%  
% %   figure()
% %   plot(rad2deg(pitch_vector), P_e3)
% %   
%   % polynomial interpolation
%   
%   [pl, S] = polyfit(pitch_vector_e3, P_e3, 3); % coefficients of the regression
%   P_fit = polyval(pl, pitch_vector_e3, S); % fit the regression
% 
%   figure()
%   plot(pitch_vector_e3, P_e3)
%   hold on
%   plot(pitch_vector_e3, P_fit)
%   yline(P_rated)
%   hold off
%   legend('Computed', 'Interpolated', 'Location','northwest')
% 
%   figure()
%   res = P_e3 - P_fit; % residuls of the regression
%   [h, p_value] = chi2gof(res);
%   scatterhist(pitch_vector_e3, res, "Direction","out","Location","SouthWest")
%   text(0, -5.5, sprintf("p-value: %f, h: %d", p_value, h), "FontSize",12);
% 
% % end
% poly_coeff = [pl(1) pl(2) pl(3) pl(4)-P_rated];
% sol = roots(poly_coeff)

disp('Question 2 done')
%% QUESTION 3
% first part of the question

% Build the velocity vector adding also the V0_rated
V0_vector_cut_in_out = linspace(V0_cutin, V0_cut_out, V0_cut_in_out_item);
V0_vector_cut_in_out(end + 1) = V0_rated; % add the rated velocity 
V0_vector_cut_in_out(end + 1) = 20; % add 20 (m/s) since it will need in Q5
V0_vector_cut_in_out = sort(V0_vector_cut_in_out); % sort the vector
V0_cut_in_out_item = V0_cut_in_out_item + 2; % increase the number of item

V0_rated_pos = find(V0_vector_cut_in_out == V0_rated); % indeces of the V0_rated used for the plot

pitch_vector_e3 = linspace(pitch_range_e3(1), pitch_range_e3(2), pitch_item_e3); % vector of pitch equally distributed in the range 

Theta_p_limit = zeros(3, V0_cut_in_out_item); % initialize the matrix to store results

for v=1:V0_cut_in_out_item % loop over differnet velocities
  P_e3 = zeros(1, pitch_item_e3);
  % chose a velocity
  V0_e3 = V0_vector_cut_in_out(v);

  % compute lambda
  omega_e3 = omega_max;
  lambda = omega_e3 * R / V0_e3;

  % chose Theta_p
  for p = 1:pitch_item_e3 %loop over different pitch
    Theta_p_e3 = pitch_vector_e3(p);

    [cp_partial, cT_partial,~,~] = cP_cT_partial(r_item_no_tip, r_vector, ...
    beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, R, lambda, ...
    Theta_p_e3, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho,0,0, i_max);

    P_e3(p) = 0.5*omega_e3*B*rho*V0_e3^2*trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial);
    
  end

  [~, max_pos] = max(P_e3);
  max_pos = max_pos;
  P_e3_part1 = P_e3(1:max_pos);
  P_e3_part2 = P_e3(max_pos:end);
  Theta_p_part1 = pitch_vector_e3(1:max_pos);
  Theta_p_part2 = pitch_vector_e3(max_pos:end);
  
  Theta_p_limit(1,v) = V0_e3; % velocity
  if V0_e3 < V0_rated + 1e-3 % velocity below the rated one
    Theta_p_limit(2, v) = 0;
    Theta_p_limit(3, v) = 0;
  else 
    Theta_p_limit(2,v) = interp1(P_e3_part1, Theta_p_part1, P_rated); % pitch angle for feathering
    Theta_p_limit(3,v) = interp1(P_e3_part2, Theta_p_part2, P_rated); % pitch angle for stall
  end

end
%%
pitch_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(Theta_p_limit(1, V0_rated_pos:end), rad2deg(Theta_p_limit(2, V0_rated_pos:end)), 'LineWidth', line_width)
hold on
plot(Theta_p_limit(1, V0_rated_pos:end), rad2deg(Theta_p_limit(3, V0_rated_pos:end)), 'LineWidth', line_width)
plot([0 V0_rated], [0 0], 'g--', 'LineWidth', line_width)
hold off
xlabel('Wind speed (m/s)')
ylabel('Pitch angle (°)')
title('Pitch angle to control the power')
legend('feathering', 'stall', 'not controlled zone', 'Location', 'southwest')
ax = gca;
ax.FontSize = font_size;
saveas(pitch_vs_V0, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\pitch_vs_V0.png','png');

% comparison of our results with the one of DTU refernce turbine pag. 33
pitch_vs_V0_reference_comparison = figure('Position', get(0, 'Screensize'));
plot(Theta_p_limit(1,:), rad2deg(Theta_p_limit(3,:)), 'LineWidth', line_width)
hold on
plot(velocity_reference, pitch_reference, 'LineWidth', line_width)
hold off
xlabel('Wind speed (m/s)')
ylabel('Pitch angle (°)')
title('Pitch angle to control the power')
legend('fathering', 'reference', 'Location','northwest')
ax = gca;
ax.FontSize = font_size;
saveas(pitch_vs_V0_reference_comparison, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\pitch_vs_V0pitch_vs_V0_reference_comparison.png','png');

%%
T_e3(p) = 0.5*B*rho*trapezoidal_integral(r_vector(1:r_item_no_tip), cT_partial);

% second part of the question 
P = size(1, V0_item); % initialize vector of power
T = size(1, V0_item); % initialize vector of thrust
cP = size(1, V0_item);
cT = size(1, V0_item);
pn_mat = zeros(V0_item, r_item); % matrix of pn (for rows different velocities, for columns r)
pt_mat = zeros(V0_item, r_item); % matrix of pt (for rows different velocities, for columns r)

for v=1:V0_cut_in_out_item % loop over differnet velocities
  V0_actual = V0_vector_cut_in_out(v);
  
  for i=1:r_item_no_tip % loop over the blade positions 
    r = r_vector(i);
    beta = beta_vector(i);
    thick = thick_vector(i);
    c = c_vector(i);
    sigma = sigma_function(c, B, r);
    
    if V0_actual < V0_rated - 1e-3 % velocities below the rated one
      lambda = lambda_opt;
      omega_actual = V0_actual * lambda / R;
      Theta_p = Theta_p_opt;

      [cp_partial, cT_partial, pt, pn] = cP_cT_partial(r_item_no_tip, r_vector, ...
        beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, R, lambda, ...
        Theta_p, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0_actual,omega_actual, i_max);
      
      cP_lower(v) = lambda*B/(R*A)*trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial);
      cT_lower(v) = B/A*trapezoidal_integral(r_vector(1:r_item_no_tip), cT_partial);
      
      P_lower(v) = cP_lower(v)*0.5*A*rho*V0_actual^3;
      T_lower(v) = cT_lower(v)*0.5*A*rho*V0_actual^2;

      V0_lower(v) = V0_actual;
      V0_upper = zeros(1, V0_cut_in_out_item - size(V0_lower, 2));

    else
      lambda = omega_max * R / V0_actual;
      Theta_p_feat = interp1(Theta_p_limit(1,:), Theta_p_limit(2,:), V0_actual);
      Theta_p_stall = interp1(Theta_p_limit(1,:), Theta_p_limit(3,:), V0_actual);
    
%       % compute the a and a_prime with the iterative method
%       [a_feat, a_prime_feat, ct_feat, cn_feat, ~] = induction_factor_convergence(a_guess, a_prime_guess, ...
%         R, r, lambda, beta, Theta_p_feat, B, sigma, aoa_mat, cl_mat, cd_mat, thick_prof, ...
%         thick, fake_zero, i_max);
%       [a_stall, a_prime_stall, ct_stall, cn_stall, ~] = induction_factor_convergence(a_guess, a_prime_guess, ...
%         R, r, lambda, beta, Theta_p_stall, B, sigma, aoa_mat, cl_mat, cd_mat, thick_prof, ...
%         thick, fake_zero, i_max);
%     
%       Vrel_stall = sqrt((V0_actual*(1 - a_stall))^2 + (omega_actual*r*(1 + a_prime_stall))^2); % (m/s)
%       Vrel_feat = sqrt((V0_actual*(1 - a_feat))^2 + (omega_actual*r*(1 + a_prime_feat))^2); % (m/s)
%   
%       pn_partial_stall(i) = 0.5*rho*c*cn*Vrel_stall^2;
%       pn_partial_feat(i) = 0.5*rho*c*cn*Vrel_feat^2;
%       pt_partial_stall(i) = r * 0.5*rho*c*ct*Vrel_stall^2;
%       pt_partial_feat(i) = r * 0.5*rho*c*ct*Vrel_feat^2;
        
      % items at velocity lower than V0_rated
      gap = size(V0_lower, 2);
      pos = v - gap; % position where to store values

      % feathering
      [cp_partial_feat, cT_partial_feat,~,~] = cP_cT_partial(r_item_no_tip, r_vector, ...
        beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, R, lambda, ...
        Theta_p_feat, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0_actual, omega_actual, i_max);
    
      cP_feat(pos) = lambda*B/(R*A)*trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial_feat);
      cT_feat(pos) = B/A*trapezoidal_integral(r_vector(1:r_item_no_tip), cT_partial_feat);

      P_feat(pos) = cP_feat(pos)*0.5*rho*A*V0_actual^3;
      T_feat(pos) = cT_feat(pos)*0.5*rho*A*V0_actual^2;

      % stall
      [cp_partial_stall, cT_partial_stall, pt, pn] = cP_cT_partial(r_item_no_tip, r_vector, ...
        beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, R, lambda, ...
        Theta_p_stall, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho, V0_actual, omega_actual, i_max);
    
      cP_stall(pos) = lambda*B/(R*A)*trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial_stall);
      cT_stall(pos) = B/A*trapezoidal_integral(r_vector(1:r_item_no_tip), cT_partial_stall);

      P_stall(pos) = cP_stall(pos)*0.5*rho*A*V0_actual^3;
      T_stall(pos) = cT_stall(pos)*0.5*rho*A*V0_actual^2;

      V0_upper(pos) = V0_actual;
    end
  end % stop looping over blade position r
  
  % store the vectors of pt and pn (which are span dependent)
  pn(end + 1) = 0; % add value correspnding to the tip
  pt(end + 1) = 0;  % add value correspnding to the tip
  pn_mat(v, :) = pn; % save the pn
  pt_mat(v, :) = pt; % save the pt
end % stop looping over different velocities


% Add the value corresponding to the rated velocity on the LH side 
[cp_partial_end, cT_partial_end, ~,~] = cP_cT_partial(r_item_no_tip, r_vector, ...
  beta_vector, thick_vector, c_vector, B, a_guess, a_prime_guess, R, lambda_opt, ...
  Theta_p_opt, aoa_mat, cl_mat, cd_mat, thick_prof, fake_zero, rho,0,0, i_max);

cP_lower(end + 1) = lambda_opt*B/(R*A)*trapezoidal_integral(r_vector(1:r_item_no_tip), cp_partial_end);
cT_lower(end + 1) = B/A*trapezoidal_integral(r_vector(1:r_item_no_tip), cT_partial_end);

P_lower(end + 1) = cP_lower(end)*0.5*A*rho*V0_rated^3;
T_lower(end + 1) = cT_lower(end)*0.5*A*rho*V0_rated^2;

% Add the value on the RH side of the rated velocity
V0_lower(end + 1) = V0_rated;

%% Plot the results

P_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_lower, P_lower, 'LineWidth', line_width);
hold on
plot(V0_upper, P_feat, 'LineWidth', line_width)
plot(V0_upper, P_stall, 'LineWidth', line_width)
hold off
xlabel('Wind velocity V0 (m/s)')
ylabel('Power (W)')
legend('Below rated velocity', 'Feathering', 'Stalling','Location', 'southeast')
title('Power output')
xlim([V0_cutin V0_cut_out])
ax = gca;
ax.FontSize = font_size;
saveas(P_vs_V0, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\P_vs_V0.png','png');


T_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_lower, T_lower, 'LineWidth', line_width);
hold on
plot(V0_upper, T_feat, 'LineWidth', line_width);
plot(V0_upper, T_stall, 'LineWidth', line_width);
hold off
xlabel('Wind velocity V0 (m/s)')
ylabel('Thrust (N)')
legend('Below rated velocity', 'Feathering', 'Stalling', 'Location','northwest')
title('Thrust force')
xlim([V0_cutin V0_cut_out])
ax = gca;
ax.FontSize = font_size;
saveas(T_vs_V0, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\T_vs_V0.png','png');

cP_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_lower, cP_lower, 'LineWidth', line_width);
hold on
plot(V0_upper, cP_feat, 'LineWidth', line_width);
plot(V0_upper, cP_stall, 'LineWidth', line_width);
hold off
xlabel('Wind velocity V0 (m/s)')
ylabel('cP')
legend('Below rated velocity', 'Feathering', 'Stalling')
title('Power coefficient cP')
xlim([V0_cutin V0_cut_out])
ax = gca;
ax.FontSize = font_size;
saveas(cP_vs_V0, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\cP_vs_V0.png','png');

cT_vs_V0 = figure('Position', get(0, 'Screensize'));
plot(V0_lower, cT_lower, 'LineWidth', line_width);
hold on
plot(V0_upper, cT_feat, 'LineWidth', line_width)
plot(V0_upper, cT_stall, 'LineWidth', line_width)
xlabel('Wind velocity V0 (m/s)')
ylabel('cT')
title('Thrust coefficient')
legend('Below rated velocity', 'Feathering', 'Stalling')
xlim([V0_cutin V0_cut_out])
ax = gca;
ax.FontSize = font_size;
saveas(cT_vs_V0, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\cT_vs_V0.png','png');


%% QUESTION 4

[pn_ashes, pt_ashes, r_ashes] = import_ashes(sensor_file_name); % load file from ashes

pt_comparison = figure('Position', get(0, 'Screensize'));
for v=1:size(V0_e4_vector, 2) % loop over different velocities
V0_e4 = V0_e4_vector(v);
pt_plt = zeros(1, r_item);
pn_plt = zeros(1, r_item);

  for i=1:r_item % loop over blade elements building the vectors to be plotted
    pt_plt(1,i) = interp1(V0_vector_cut_in_out, pt_mat(:,i), V0_e4);
  end

% do subplots for different velocities
subplot(2,2,v)
plot(r_vector, pt_plt);
hold on 
plot(r_ashes, pt_ashes(v, :))
xlabel('r (m)')
ylabel('pt')
legend('BEM code', 'Ashes')
title(strcat('V_0 = ', num2str(V0_e4), ' (m/s)'))
set(gca, 'FontSize', font_size/2);
end
sgtitle('p_t coefficients')

saveas(pt_comparison, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\pt_comparison.png','png');


pn_comparison = figure('Position', get(0, 'Screensize'));
for v=1:size(V0_e4_vector, 2) % loop over different velocities
V0_e4 = V0_e4_vector(v);
pt_plt = zeros(1, r_item);
pn_plt = zeros(1, r_item);

  for i=1:r_item % loop over blade elements building the vectors to be plotted
    pn_plt(1,i) = interp1(V0_vector_cut_in_out, pn_mat(:,i), V0_e4);
  end

% do subplots for different velocities
subplot(2,2,v)
plot(r_vector, pn_plt)
hold on 
plot(r_ashes, pn_ashes(v, :))
xlabel('r (m)')
ylabel('pn')
legend('BEM code', 'Ashes')
title(strcat('V_0 = ', num2str(V0_e4), ' (m/s)'))
set(gca, 'FontSize', font_size/2);
end
sgtitle('p_n coefficients')
saveas(pn_comparison, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\pt_comparison.png','png');




%% QUESTION 5
% first part
P_e5_feat = [P_lower(1:end-1), P_feat]; % power vector for velocities from cut in up to cut out
P_e5_stall = [P_lower(1:end-1), P_stall]; % power vector for velocities from cut in up to cut out

V0_e5 = [V0_lower(1:end-1), V0_upper]; % velocity vector 

AEO = 0;  % Initialize AEO
AEO_20 = 0; % 

[AEO_25_feat, AEO_20_feat, CF_25_feat, CF_20_feat, delta_AOE_feat] = Annual_Energy_Output(P_e5_feat, V0_e5, P_rated, A_weibull, k_weibull);

[AEO_25_stall, AEO_20_stall, CF_25_stall, CF_20_stall, delta_AOE_stall] = Annual_Energy_Output(P_e5_stall, V0_e5, P_rated, A_weibull, k_weibull);

display(strcat('When feathering the energy difference is ', num2str(delta_AOE_feat/1e9), ' (GWh)'));
display(strcat('When stalling the energy difference is ', num2str(delta_AOE_stall/1e9), ' (GWh)'));

V0_weibull = [0:1:V0_cut_out];
weibull = zeros(1, size(V0_weibull, 2));
for i=1:size(V0_weibull, 2)
  weibull(i) = weibull_pdf(V0_weibull(i), A_weibull, k_weibull);
end

weibull_plot = figure('Position', get(0, 'Screensize'));
plot(V0_weibull, weibull, 'LineWidth', line_width) 
xlabel('Wind speed (m/s)')
ylabel('Weibull PDF')
title('Weibull pdf as function of the wind velocity')
ax = gca;
ax.FontSize = font_size;
saveas(weibull_plot, 'C:\Users\Niccolò\Documents\UNIVERSITA\5° ANNO\WIND_ENERGY\assignment_1\figures\weibull_plot.png','png');