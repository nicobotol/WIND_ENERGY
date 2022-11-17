function [CL,CD,CM] = ...
    extract_from_matrix(aerodynamics, positions)

    %% Description
    % This function extracts the lift, drag and momentum coefficients from
    % the matrix containing them
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
    % - aerodynamics  ---> Table with all the coefficients
    % - positions     ---> Positions of the coefficients in the table
    % _____________________________________________________________________
    % OUTPUTS
    % - CL ---> Array with the lift coefficient
    % - CD ---> Array with the lift coefficient
    % - CM ---> Array with the lift coefficient
    % _____________________________________________________________________
    % AUXILIARY
    % - cl_col        ---> Position of cl in the table
    % - cd_col        ---> Position of cd in the table
    % - cm_col        ---> Position of cm in the table
    % - N_alpha       ---> Number of angles of attack
    % - N_data        ---> Number of different coefficients
    % - N_tc          ---> Number of thicknesses (different)
    % - ii            ---> Iteration variable 1 for "for loops" [-]
    % _____________________________________________________________________
    % *********************************************************************
    %% Operation
    % *********************************************************************

    N_data = size(positions,2);
    cl_col = positions(1);
    cd_col = positions(2);
    cm_col = positions(3);

    N_alpha = size(aerodynamics, 1);       % Number of aoa
    N_tc  = size(aerodynamics,2) / N_data; % Number of thicknesses
    
    CL = zeros(N_alpha, N_tc);
    CD = zeros(N_alpha, N_tc);
    CM = zeros(N_alpha, N_tc);
    
    for ii = 1:N_tc
        CL(:,ii) = aerodynamics(:,cl_col + N_data*(ii-1));
        CD(:,ii) = aerodynamics(:,cd_col + N_data*(ii-1));
        CM(:,ii) = aerodynamics(:,cm_col + N_data*(ii-1));
    end

end