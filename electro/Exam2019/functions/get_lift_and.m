function [Cl, Cd, Cm] = get_lift_and(CL, CD, CM, aoa_ser, aoa_targ)
    %% Description
    % This function interpolates the lift, drag and momentum coefficients
    % for the target angle of attack
    %% General Information
    % Version: 1
    % Date: 20/09/2021
    % Authors: Sowmya, Philipp, Carlos
    % Denmark Technical University (DTU)
    % Wind Turbine Technologies and Aerodynamics
    % Assignment 1
    %% Function dictionary
    % _____________________________________________________________________
    % INPUTS
    % - CL       ---> Array with the lift coefficient for each thickness [-]
    % - CD       ---> Array with the drag coefficient "  " " [-]
    % - CM       ---> Array with the moment coefficient " " " [-]
    % - aoa_ser  ---> Series of angle of attack for CL,CD,CM [deg]
    % - aoa_targ ---> Target angle of attack for interpolation [deg]
    % _____________________________________________________________________
    % OUTPUTS
    % - Cl ---> Lift coefficient for the desired tc and alpha [-]
    % - Cd ---> Drag coefficient for the desired tc and alpha [-]
    % - Cm ---> Moment coefficient for the desired tc and alpha [-]
    % _____________________________________________________________________
    % AUXILIARY
    % - A             ---> Auxiliary matrix for storing the data [-]
    % - AA            ---> Auxiliary array for comparing values [-]
    % - alpha         ---> Angles of attack of the aerodyn. data [deg]
    % - blade_tc      ---> Vector with the *different* tc values of blade
    % - cd            ---> Drag coefficient of the aerodyn. data [-]
    % - cl            ---> Lift coefficient of the aerodyn. data [-]
    % - cm            ---> Momentum coefficient of the aerodyn. data [-]
    % - cl_int        ---> Lift coefficient interpolated [-]
    % - cd_int        ---> Drag coefficient interpolated [-]
    % - cm_int        ---> Momentum coefficient interpolated [-]
    % - eps           ---> Small value for comparison [-]
    % - ii            ---> Iteration variable 1 for loops [-]
    % - kk            ---> Iteration variable 2 for loops [-]
    % - N_aoainterp   ---> Number of aoa for interpolation [-]
    % - N_cols        ---> Number of columns per thickness (cl, cd, cm) [-]
    % - N_files       ---> Number of files of airfoils provided [-]
    % - N_tc          ---> Number of different tc ratios in the blade [-]
    % - not_found     ---> Auxiliary boolena variable [bool]
    % _____________________________________________________________________
    % *********************************************************************
    %% Operation
    % *********************************************************************
    
    Cl = interp1(aoa_ser, CL, aoa_targ);
    Cd = interp1(aoa_ser, CD, aoa_targ);
    Cm = interp1(aoa_ser, CM, aoa_targ);
end