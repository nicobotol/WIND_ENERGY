function [radius, chord, twist, thickness] = input_bem_data(filename)
    %% Description
    % This function reads the data from a file with the following structure
    % | radius | twist angle | chord | thickness |
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
    % - Filename  ---> String with the absolute or relative path to the
    %                   file
    % _____________________________________________________________________
    % OUTPUTS
    % - radius    ---> Radial position of the blade element
    % - chord     ---> Airfoil chord of the blade element
    % - twist     ---> Twist angle of the blade element
    % - thickness ---> Thickness ratio of the blade element
    % _____________________________________________________________________
    % AUXILIARY
    % - A         ---> Auxiliary matrix for storing the data
    % _____________________________________________________________________
    % *********************************************************************
    %% Operations
    % *********************************************************************
    
    % Import the data
    A = readmatrix(filename);
    
    % Save the data in output arrays
    radius    = A(:,1);                 % Radial position of the elements
    twist     = A(:,2);                 % Twist angle
    chord     = A(:,3);                 % Chord
    thickness = A(:,4);                 % Thickness ratio
end