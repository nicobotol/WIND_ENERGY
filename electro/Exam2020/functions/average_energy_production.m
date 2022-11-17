function AEP = average_energy_production(Power, V,k,A)
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
    % - Power   ---> Values of the power curve [W]
    % - V       ---> Series of wind speed [m/s]
    % - k       ---> Shape factor [-]
    % - A       ---> Scaling factor [-]
    %
    % _____________________________________________________________________
    % OUTPUTS
    % - f       ---> Probability distribution [-]
    % _____________________________________________________________________
    % AUXILIARY
    % - N_Power ---> Number of power points [-]
    % - AEP_ii  ---> Auxiliary variable to accumulate the AEP values [Wh]
    % _____________________________________________________________________

    % Initialize the Average Energy Production for integration
    AEP = 0;
    N_Power = length(Power);
    for ii = 1:N_Power-1
        AEP_ii      = 0.5 * (Power(ii) + Power(ii+1)) ... % Power
                        * 8760 ...                      % Hours of the year
                        * weibull(V(ii),k,A);           % Probability
                    
        AEP         = AEP + AEP_ii;                     % Cumulative
    end
    
end