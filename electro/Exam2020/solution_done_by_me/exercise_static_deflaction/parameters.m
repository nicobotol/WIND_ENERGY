%% PARAMETERS

font_size = 32; % font szie for plots
line_width = 2; % line width for plots

i_max = 300;  % maximum number of accepted iterations
fake_zero = 1e-8; % thershold for exiting the loop

V0 = 11; % windspeed to investigate
R = 89.17; % rotor radius (m)
A = R^2 * pi; % rotor area (m^2)
B = 3; % number of blades (#)
rho = 1.225;  % air density (kg/m^3)
V0_cutin = 4; % cut in wind velocity (m/s)
V0_cut_out = 25; % cut out wind velocity (m/s)
V0_rated = 11.44; % rated wind speed
P_rated = 10.64e6; % rated power (W)
omega_max = 1.01; % rated rotational speed (rad/s)
cp_opt = 0.465; % optimal cp
%lambda_opt = 7.857; % optimal lambda
lambda_opt = 7.5;
Theta_p_opt = 0*pi/180; % optimla pitch angle [rad]

lambda_range = [5 10]; % range in within look for optimal lambda
pitch_range = [deg2rad(-2) deg2rad(5)]; % range in within look for optimal picth angle
lambda_item = 15; % number of divisions of lambda range (initially it was 10)
pitch_item = 15 ; % number of divisions of pitch range (initially it was 10)

V0_item = 25; % number of divison for the omega plot in Q2
V0_cut_in_out_item = 30; % number of division for velocity vector in Q3
pitch_item_e3 = 30; % number of division for the pitch in Q3
pitch_range_e3 = [-0.22 0.42];  % (rad) range in within look for picth angle 

V0_e4_vector = [5, 9, 11, 20];

% parameters of weibull distribution 
A_weibull = 9;
k_weibull = 1.9;

% first guess for a and a_prime
a_guess = 0;
a_prime_guess = 0.1;

% file with airfoils parameter 
filenames = [ "airfoil_data\cylinder", "airfoil_data\FFA-W3-600",  "airfoil_data\FFA-W3-480", "airfoil_data\FFA-W3-360", "airfoil_data\FFA-W3-301", "airfoil_data\FFA-W3-241"];
% file with blade parameters
blade_filename = "airfoil_data\bladedat.txt";
% t/c ratio (they are in the same order provided as the file uploaded)
thick_prof = [100 60 48 36 30.1 24.1];

% pitch parameters from DTU reference turbine pag 33
velocity_reference = [4:1:25];
pitch_reference = [2.751 1.966 0.896 0.000 0.000 0.000 0.000 0.000 4.502 7.266 9.292 10.958 12.499 13.896 15.200 16.432 17.618  18.758 19.860 20.927 21.963 22.975 ];

% Ashes comparison
V0_ashes = [5, 9, 11, 20]; % velocities
P_ashes = [0.894, 5.14, 9.4, 10.8]; % power(MW)
T_ashes = [311, 1012, 1515, 668]; % thrust (kN)

