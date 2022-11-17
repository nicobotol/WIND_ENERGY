function idx_tc = which_tc(tc_vec, tc_be)
%% Description
    % This function returns the value of the thickness
    %% General Information
    % Version: 1
    % Date: 19/09/2021
    % Authors: Sowmya, Philipp, Carlos
    % Denmark Technical University (DTU)
    % Wind Turbine Technologies and Aerodynamics
    % Assignment 1
    %% Function dictionary
    % _____________________________________________________________________
    % INPUTS
    % - tc_vec ---> Array with the different thicknesses of the blade [-]
    % - tc_be  ---> Thickness of the target blade element [%]
    % _____________________________________________________________________
    % OUTPUTS
    % - idx_tcc ---> Position of the thickness in the tc_vec array
    % _____________________________________________________________________
    % AUXILIARY
    % - AA      ---> Auxiliary array for calculating minimum value[-]
    
    
    AA = tc_vec - tc_be;
    [~,idx_tc] = min(abs(AA));
    
end