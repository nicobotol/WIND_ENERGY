function pitch_control = ...
    calculate_pitch(c_p, pitch_vals, c_p_target, P_R, rho, R, V0)
%% Description
    % This function calculates the values of pitch required to keep the
    % power below rated for high wind speeds
    %% General Information
    % Version: 1
    % Date: 26/09/2021
    % Authors: Sowmya, Philipp, Carlos
    % Denmark Technical University (DTU)
    % Wind Turbine Technologies and Aerodynamics
    % Assignment 1
    %% Function dictionary
    % _____________________________________________________________________
    % INPUTS
    % - c_p ---> Matrix with C_P values from BEM [-]
    % - c_T ---> Matrix with C_T values from BEM [-]
    % - lambda_vals  ---> TSR values for the C_P Matrix [-]
    % - pitch_vals  ---> Theta (pitch) values for the C_P Matrix [-]
    % - lambda  ---> TSR target [-]
    % - V0  ---> Wind speed target [m/s]
    % _____________________________________________________________________
    % OUTPUTS
    % - 
    % _____________________________________________________________________
    % AUXILIARY
    % - 
    

    % Calculate the maximum power coefficient to keep power below rated
    c_pmax = P_R / (0.5 * rho * V0^3 * pi * R^2);   
    
    N_Theta = size((pitch_vals),2);
    
    [~, idx_cppeak] = max(c_p);
    cp_lhs = zeros(idx_cppeak,1);
    cp_rhs = zeros(N_Theta - idx_cppeak + 1,1);
    pitch_vals_lhs = zeros(idx_cppeak,1);
    pitch_vals_rhs = zeros(N_Theta - idx_cppeak + 1,1);
    if idx_cppeak > 1
        cp_lhs(:) = c_p(1:idx_cppeak);
        cp_rhs(:) = c_p(idx_cppeak:N_Theta);
        pitch_vals_lhs(:) = pitch_vals(1:idx_cppeak);
        pitch_vals_rhs(:) = pitch_vals(idx_cppeak:N_Theta);
        theta_p1 = find_theta(cp_lhs, pitch_vals_lhs, c_pmax);
    else
        cp_rhs = c_p_target;
        pitch_vals_rhs = pitch_vals(idx_cppeak:N_Theta);
        theta_p1 = -99;
    end
    
    theta_p2 = find_theta(cp_rhs, pitch_vals_rhs, c_pmax);
    pitch_control = [theta_p1 theta_p2];
    
end