function vec_different = get_different(values, tolerance)

    %% Description
    % This function gives the different values encountered in a vector with
    % a certain tolerance
    %% General Information
    % Authors: Sowmya, Philipp, Carlos
    % Denmark Technical University (DTU)
    % Wind Turbine Technologies and Aerodynamics
    % Assignment 1
    %% VERSIONS
    % _____________________________________________________________________
    % Version: 1 Date: 19/09/2021
    % Write the function and first tests
    % _____________________________________________________________________
    % Version: 2 Date: 04/10/2021
    % Sort the values before searching for different, for cases where the
    % blade thickness is not purely descending
    % 
    %% Function dictionary
    % _____________________________________________________________________
    % INPUTS
    % - values    ---> 1D array
    % - tolerance ---> Tolerance for the comparison
    % _____________________________________________________________________
    % OUTPUTS
    % - vec_different ---> Vector with the different values
    % _____________________________________________________________________
    % AUXILIARY
    % - aux_1         ---> Auxiliary scalar 1 for comparisons
    % - aux_2         ---> Auxiliary scalar 2 for comparisons
    % - count         ---> Auxiliary variable for counting
    % - ii            ---> Iteration variable 1 for "for loops" [-]
    % _____________________________________________________________________
    % *********************************************************************
    %% Operation
    % *********************************************************************
    % Count different thicknesses in the main file
    
    N_vals = size(values,1);
    
    aux_2 = -99; % Thickness we will never have (for first comparison
    count = 0;
    values = sort(values, 'descend');
    for ii = 1:N_vals
        aux_1 = values(ii);
        if abs(aux_1 - aux_2) < tolerance
            continue
        else
            aux_2 = aux_1;
            count = count+1;
        end
    end
    % Pass the quantity of different values to a new variable
    number_different = count;
    % Preallocate the output vector
    vec_different = zeros(number_different,1);
    % Restart the count and the first value for comparison
    count = 0;
    aux_2 = -99;
    for ii = 1:N_vals
        aux_1 = values(ii);
        if abs(aux_1 - aux_2) < tolerance
            continue
        else
            aux_2 = aux_1;
            count = count+1;
            vec_different(count) = aux_2;
        end
        
    end
end