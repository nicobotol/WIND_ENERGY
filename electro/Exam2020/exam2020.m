%% Exam 2020

%% Q1 · VAWT

Wx = 4; Wy = 0; V0 = 12; omega = 115*pi/30; R = 2;
theta = [0 pi/2 pi 3*pi/2];

% Coordinate Shifting
x           = -R*sin(theta);
y           = R*cos(theta);

% Velocities
Vrel_x      = omega*y + V0 - Wx;
Vrel_y      = -omega*x + Wy;
Vnorm       = (V0-Wx)*sin(theta) - Wy*cos(theta);
Vtan        = (V0-Wx)*cos(theta) + Wy*sin(theta) + omega*R;
Vrel        = sqrt(Vnorm.^2 + Vtan.^2);

% Angle of attack
aoa         = atand(Vnorm./Vtan);

%% Q2 · 
T = 800e3; V0 = 10; rho = 1.225; R = 70;
CT = T / (0.5 * rho * pi * R^2 * V0^2);

syms a;
eq1 = CT - 4*a*(1-a) == 0;
sol1 = double(solve(eq1,a));
eq  = CT - 4*a*(1 - 0.25*(5-3*a)*a) == 0;
sol = double(solve(eq,a));

a_w = min(sol1);

Cp = 4*a_w*(1-a_w)^2;
P = Cp * 0.5 * rho * pi * R^2 * V0^3;

%% Q3 · Pitch for rated power

V0 = 15; theta_p = (0:0.5:25)'; 
omega_r = 1.2; B = 3; R  = 89.17; glc = 1;
rho_air = 1.225; 

functions_folder = ".\functions"; input_folder = ".\inputs";
addpath(functions_folder);  addpath(input_folder);
file_aoa_ser    = "aoa_series.dat"; aoa_ser    = importdata(file_aoa_ser);
file_bem_data   = "bladedat.txt"; 
[r, c, beta, tc] = input_bem_data(file_bem_data); NE = size(r,1);
file_blade_aero = "aerodynamics_bem.dat"; positions = [1 2 3];
file_performance = "BEM_performance.mat";
% The aerodynamic files should be in increasing distance to hub
files_aero = ["cylinder.txt" "FFA-W3-600.txt" "FFA-W3-480.txt" ...
    "FFA-W3-360.txt" "FFA-W3-301.txt" "FFA-W3-241.txt"];
tc_files     = [100, 60, 48, 36, 30.1, 24.1]; % Thicknesses ratios
% Import aerodynamic data
table_aerodynamics = importdata(file_blade_aero);
[CL,CD,CM] = extract_from_matrix(table_aerodynamics, positions);

lambda = omega_r * R / V0;
cp  = zeros(length(theta_p),1);
PWT = zeros(length(theta_p),1);

for ii = 1:length(theta_p)
    % Call BEM function
    [cp(ii), ~,~,~,~,~,~,~] = ...
        BIG_BEMalg(R,B,rho_air,V0,lambda,theta_p(ii), CL, CD, CM, aoa_ser, ...
        glc, r, c, beta, tc);

    PWT(ii) = 0.5*cp(ii)*rho_air*pi*R^2*V0^3;
end

theta_target = interp1(PWT, theta_p, 10e+6);
figure('Name', 'find pitch')
plot(theta_p, PWT)

%% Q4 · Tangential induction factor
clear all; close all; clc;

u  = 6; theta_p = deg2rad(1); r = 15; R = 30; lambda = 7;
V0 = 9; alpha = deg2rad(2); beta = deg2rad(7);
omega = lambda * V0 / R;
phi = alpha + beta + theta_p;
a = 1 - u/V0;
a_prime = 1/(omega*r*tan(phi)) * (1 - a)*V0 - 1;

%% Q5 · 
clear all; close all; clc;
F = load('flexibility.mat', 'F');
F = F.F;
N_F = size(F,1);
loads = zeros(N_F,1);
loads(1:2:end-1) = 0;

pz = 0;
found = 0; counter = 0; 
while ~found

    loads(2:2:end) = pz;
    u = F * loads;
    if u(end) > 5
        break
    else
        pz = pz + 1;
        counter = counter + 1;
        p_z(counter) = pz;
        u_target(counter) = u(end);
    end
end

%% Q8 · 
A = 8; k = 2;
V = 5:0.01:25;

for ii = 1:length(V)
    V0 = V(ii);
    if V0 < 10
        P(ii) = 0.7e6;
    elseif V0 < 13
        P(ii) = 2e6;
    else
        P(ii) = 3e6;
    end
end

AEP = average_energy_production(P, V, k, A);
