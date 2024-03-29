%% PARAMETERS

i_max = 300;  % maximum number of accepted iterations
fake_zero = 1e-8; % thershold for exiting the loop

R = 89.17; % rotor radius (m)
A = R^2 * pi; % rotor area (m^2)
B = 3; % number of blades (#)
rho = 1.225;  % air density (kg/m^3)
V0_cutin = 4; % cut in wind velocity (m/s)
V0_cut_out = 25; % cut out wind velocity (m/s)
P_rated = 10.64e6; % rated power (W)

lambda_range = [5 10]; % range in within look for optimal lambda
pitch_range = [deg2rad(-2) deg2rad(5)]; % range in within look for optimal picth angle
lambda_item = 10; % number of divisions of lambda range
pitch_item = 10; % number of divisions of pitch range

V0_item = 20; % number of divison for the omega plot in Q2
V0_cut_in_out_item = 30; % number of division for velocity vector in Q3

%V0 = 8;
%omega = 2.61;
%Theta_p = -3*pi/180;
%beta = 2*pi/180;
%c = 1.5;
%cl = 0.5;
%cd = 0.01;
%r = 24.5;

% first guess for a and a_prime
a_guess = 0;
a_prime_guess = 0.1;

% file with airfoils parameter 
filenames = [ "airfoil_data\cylinder", "airfoil_data\FFA-W3-600",  "airfoil_data\FFA-W3-480", "airfoil_data\FFA-W3-360", "airfoil_data\FFA-W3-301", "airfoil_data\FFA-W3-241"];
% file with blade parameters
blade_filename = "airfoil_data\bladedat.txt";
% t/c ratio (they are in the same order provided as the file uploaded)
thick_prof = [100 60 48 36 30.1 24.1];

