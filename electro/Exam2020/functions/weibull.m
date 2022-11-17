function f = weibull(v,k,A)
    %% Description
    % This function calculates the Weibull distribution given the wanted
    % wind speed vector, the shape factor and the scaling factor
    %% General Information
    % Version: 1
    % Date: 29/09/2021
    % Authors: Sowmya, Philipp, Carlos
    % Denmark Technical University (DTU)
    % Wind Turbine Technologies and Aerodynamics
    % Assignment 1
    %% Function dictionary
    % _____________________________________________________________________
    % INPUTS
    % - v       ---> Series of wind speed [m/s]
    % - k       ---> Shape factor [-]
    % - A       ---> Scaling factor [-]
    %
    % _____________________________________________________________________
    % OUTPUTS
    % - f       ---> Probability distribution [-]
    % _____________________________________________________________________

    
    f = k / A .* (v / A) .^(k-1) .* exp(-(v/A).^k);
    
end