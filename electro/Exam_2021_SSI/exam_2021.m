%% Wind Turbine Technology Exam 2021
% Author: Sowmya Srinivasan Iyer
% Student number: S213852

clear all; close all; clc;

%% BEM Question
clear all; close all; clc;

V0 = 10; theta_p = 0; omega_r = 0.9; B = 3; R  = 25; glc = 1;
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

% Call BEM function
[cp, cT,MB_normal,MB_tangential, p_normal, p_tangential,...
    T_normal, T_tangential] = ...
    BIG_BEMalg(R,B,rho_air,V0,lambda,theta_p, CL, CD, CM, aoa_ser, ...
    glc, r, c, beta, tc);

PWT = 0.5*cp*rho_air*pi*R^2*V0^3;

% Plot the loads
figure('Name','Loads')
plot(r, p_normal); hold on;
plot(r, p_tangential); grid on;

% Q6 · Iterate for values of omega

c_new = c .* 0.7;

omega_vals = (0:0.01:2)';
NOmega_vals = length(omega_vals);
P = zeros(NOmega_vals, 1);
for ii = 1:NOmega_vals
    omega_r = omega_vals(ii);
    lambda = omega_r * R / V0;
    [cp, cT,MB_normal,MB_tangential, p_normal, p_tangential, ...
            T_normal, T_tangential] = ...
        BIG_BEMalg(R,B,rho_air,V0,lambda,theta_p, CL, CD, CM, aoa_ser, ...
        glc,r,c_new,beta,tc);
    P(ii) = cp * 0.5 * rho_air * pi * R^2 * V0^3;
end

figure('Name','Power curve')
plot(omega_vals, P); hold on;
yline(PWT); grid on

%% Question 1: 
clear all; close all; clc;
V0 = 9; %m/s
R = 2; %m
CT = 0.7;
rho = 1.225; %kg/m^3

T = 0.5*CT*rho*V0^2*pi*R^2;


syms a;
eq1 = CT - 4*a*(1-a) == 0;
sol1 = double(solve(eq1,a));
eq  = CT - 4*a*(1 - 0.25*(5-3*a)*a) == 0;
sol = double(solve(eq,a));

a_w = min(sol1);
a_final = 0.2261

u = 1 + a_final*V0;





%% Question 2: 
clear all; close all; clc;
V0 = 9; omega = (0:1:5000)*pi/30; R = 2;
theta = pi;

% Coordinate Shifting
x           = -R;
y           = 0;

% Velocities

% Vnorm       = (V0-Wx)*sin(theta) - Wy*cos(theta);
% Vtan        = (V0-Wx)*cos(theta) + Wy*sin(theta) + omega*R;
% Vrel        = sqrt(Vnorm.^2 + Vtan.^2);
Vrel = 6;

% Angle of attack
aoa         = deg2rad(6);

for ii = 1:length(omega)
    Vrel_x(ii)      = omega(ii)*y + V0 - V0*cos(aoa);
    Vrel_y(ii)      = -omega(ii)*x + V0*sin(aoa);
    Vrel_test(ii)       = sqrt(Vrel_x(ii)^2 + Vrel_y(ii)^2);
end


omega_target = interp1(Vrel_test,omega,Vrel);

%% Question 3: 
clear all; close all; clc;

R = 25; B = 2;
r = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 25 25];
pT = [0 0 0 0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 136 102 68 24 0];
Moment = r.^2 .*pT;
torque = B*sum(Moment);




% 
% clear all; close all; clc;
% 
% V0 = 15; theta_p = 10; 
% omega_r = 12.2; B = 2; R  = 25; glc = 1;
% rho_air = 1.225; 
% 
% functions_folder = ".\functions"; input_folder = ".\inputs";
% addpath(functions_folder);  addpath(input_folder);
% file_aoa_ser    = "aoa_series.dat"; aoa_ser    = importdata(file_aoa_ser);
% file_bem_data   = "bladedat.txt"; 
% [r, c, beta, tc] = input_bem_data(file_bem_data); NE = size(r,1);
% file_blade_aero = "aerodynamics_bem.dat"; positions = [1 2 3];
% file_performance = "BEM_performance.mat";
% % The aerodynamic files should be in increasing distance to hub
% files_aero = ["cylinder.txt" "FFA-W3-600.txt" "FFA-W3-480.txt" ...
%     "FFA-W3-360.txt" "FFA-W3-301.txt" "FFA-W3-241.txt"];
% tc_files     = [100, 60, 48, 36, 30.1, 24.1]; % Thicknesses ratios
% % Import aerodynamic data
% table_aerodynamics = importdata(file_blade_aero);
% [CL,CD,CM] = extract_from_matrix(table_aerodynamics, positions);
% 
% 
% lambda = omega_r * R / V0;
% cp  = zeros(length(omega_r),1);
% PWT = zeros(length(omega_r),1);
% 
% for ii = 1:length(omega_r)
%     % Call BEM function
%     [cp(ii), cT,MB_normal,MB_tangential, p_normal, p_tangential(ii), ...
%             T_normal, T_tangential] = ...
%         BIG_BEMalg(R,B,rho_air,V0,lambda,theta_p, CL, CD, CM, aoa_ser, ...
%         glc, r(ii), c, beta, tc);
% 
%     PWT(ii) = 0.5*cp(ii)*rho_air*pi*R^2*V0^3;
%     Torque(ii) = p_tangential(ii)*r(ii);
% end
% 
% omega_target = interp1(PWT, omega, 10e+6);
% figure('Name', 'find omega')
% plot(omega_r, PWT)

%% Question 4: 
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
%% Question 5: 
clear all; close all; clc;

%% Question 6:
clear all; close all; clc;

const = 1:1:100;
for MM = 1:length(const)
    % Static deflection and modeshapes
    count = 1;
    for ii = 1:2:2*(N)
        m(ii,ii) = const*mass(count);
        m(ii+1,ii+1) = const*mass(count);

        count = count+1;
    end

    F = zeros(2*N,2*N);    % Flexibility Matrix


    for ii = 1:2*N

        p       = zeros(2*N,1);
        p(ii)   = 1;
        py      = p(1:2:end);
        %     py = 0;
        pz      = p(2:2:end);
        %     pz = 4816;

        statics_OUT = statics_fun(bladedata,py,pz,0);

        uy = statics_OUT.uy;
        uz = statics_OUT.uz;

        u            = zeros(2*N,1);
        u(1:2:end)   = uy;
        u(2:2:end)   = uz;

        F(:,ii)      = u;
    end


    [V,D] = eigs(F*m);
    omega_1(MM) = sqrt(1/D(1,1));
    omega_2 = sqrt(1/D(2,2));
    omega_3 = sqrt(1/D(3,3));

end

%% Question 7: 
clear all; close all; clc;

