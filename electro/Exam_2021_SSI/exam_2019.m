%%% 2019 EXAM
clear all; close all; clc;
%% Q1 · Calculate thrust coefficient
% Parameters
V0 = 9; T = 1000e3; rho_air = 1.225; R = 70;
% Calculate thrust coefficient
CT = T / (0.5 * rho_air * pi * R^2 * V0^2);

%% Q2 · Estimate axial induction factor for same conditions
syms a;
eq = CT - 4*a*(1 - 0.25*(5-3*a)*a) == 0;
sol = double(solve(eq,a));

%% Q3 · Formula to compute the power from loads in HAWT
% P = omega * M_R = omega * B * integral of (r * pt)

%% Q4 · VAWT
R = 2; % Radius in meters
omega_R = 115*pi/30; %Speed in rpm
B = 3; % Number of blades
pitch = 0; % Pitch angle in degrees
c = 0.1; % Chord in meters
V0 = 7; dt = 0.001; ntime = 1e4; target_angles = [-pi/2 0 pi/2 pi];
% Call VAWT solver
[t,px_tot,py_tot,pt_tot,CP,CT, a, theta, wx, wy] = ...
                        VAWT(B, dt, ntime, R, V0, omega_R, c, rho_air);

% Let's do it the easy way
Wg = a(end)*V0; Wx = Wg*(1-0.3*sin(target_angles));
Wy = 0.3*Wx.*cos(target_angles);
% Truncate data for the last values
aux = theta(end,3); aux2 = aux - 2*pi;
[last_min_theta, idx_min_theta] = min(abs(theta(:,3) - aux2));
if last_min_theta > aux2
    idx_min_theta = idx_min_theta - 1;
    last_min_theta = theta(idx_min_theta, 3);
end
Theta = mod(theta(idx_min_theta:end, 3),2*pi); 
[Theta_max, idx_max_theta] = max(Theta);
Theta(1:idx_max_theta) = Theta(1:idx_max_theta) - 2*pi;
Wx = wx(idx_min_theta:end,3); 
figure('Name', 'Induced wind in x direction'); plot(Theta, Wx); grid on;
Wx_target = interp1(Theta, Wx, target_angles);

%% Q5 · Calculate aerodynamic loads
clear all; close all; clc;

V0 = 10; theta_p = 0; omega_r = 0.9; B = 3; R  = 89.17; glc = 1;
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
    [cp, ~, ~, ~, ~, ~, ~, ~] = ...
        BIG_BEMalg(R,B,rho_air,V0,lambda,theta_p, CL, CD, CM, aoa_ser, ...
        glc,r,c_new,beta,tc);
    P(ii) = cp * 0.5 * rho_air * pi * R^2 * V0^3;
end

figure('Name','Power curve')
plot(omega_vals, P); hold on;
yline(PWT); grid on

%% Q7 · Bending moment 


%% Functions

function [t,px_tot,py_tot,pt_tot,CP,CT, a, theta, Wx, Wy] = ...
                                    VAWT(B,dt,ntime, R, V0, omega, c, rho)

    % Data Input
    file_name = 'airfoil.txt';
    airfoil   = importdata(file_name);

    aoa         = deg2rad(airfoil.data(:,1));           % AoA given in DEG, then to RAD
    lift_coeff  = airfoil.data(:,2);
    drag_coeff  = airfoil.data(:,3);

    % Parameters
    tau     = 2*R/V0;   % Time Relaxation constant   [s]

    a       = zeros(ntime,1);
    a(1)    = 0; % Initialize a=0. Redundant but just for protocol.

    for n = 1:ntime
        
        t(n)        = (n-1)*dt;
        theta1      = omega*t(n);
        W_global    = a(n) * V0;
        
        for ii = 1:B    % for each blade
            theta(n,ii) = theta1 + 2*pi*(ii-1)/B; % Shifting Theta per Blade
            azi(n,ii)   = mod(theta(n,ii),2*pi);
            Wx(n,ii)    = W_global * (1 - 0.3*sin(theta(n,ii)));
            Wy(n,ii)    = 0.3*Wx(n,ii)*cos(theta(n,ii));
            
            % Coordinate Shifting
            x           = -R*sin(theta(n,ii));
            y           = R*cos(theta(n,ii));
            
            % Velocities
            Vrel_x      = omega*y + V0 - Wx(n,ii);
            Vrel_y      = -omega*x + Wy(n,ii);
            Vnorm       = (V0-Wx(n,ii))*sin(theta(n,ii)) - Wy(n,ii)*cos(theta(n,ii));
            Vtan        = (V0-Wx(n,ii))*cos(theta(n,ii)) + Wy(n,ii)*sin(theta(n,ii)) + omega*R;
            Vrel        = sqrt(Vnorm^2 + Vtan^2);
            
            alpha(n,ii) = atan(Vnorm/Vtan);
            pitch       = 0;                                % NO PITCH!
            phi         = alpha(n,ii) + pitch;
            
            % Interpolate
            Cl(n,ii)    = interp1(aoa,lift_coeff,alpha(n,ii));
            Cd(n,ii)    = interp1(aoa,drag_coeff,alpha(n,ii));
            
%             c           = S*R/B;
            
            % Forces
            l           = 0.5*rho*Vrel^2*c*Cl(n,ii);
            d           = 0.5*rho*Vrel^2*c*Cd(n,ii);
            
            % Loads
            cos_beta    = Vrel_y / Vrel;
            sin_beta    = Vrel_x / Vrel;
            
            p_x(n,ii)   = l * cos_beta + d * sin_beta;
            p_y(n,ii)   = -l* sin_beta + d * cos_beta;
            
            p_n(n,ii)   = l*cos(phi) + d*sin(phi);
            p_t(n,ii)   = l*sin(phi) - d*cos(phi);
        end
        
        % Thrust and Power Coefficients
        % Total
        CT(n) = sum(p_x(n,:)) / (rho*V0^2*R);
        CP(n) = omega*sum(p_t(n,:)) / (rho*V0^3);
        
        
        % Vary a
        if a(n) <= 1/3
            fg = 1;
        else
            fg = 1/4 * (5-3*a(n));                  % Quasi-Steady induced Wind
        end
        
        a_qs   = CT(n) / (4*(1-fg*a(n)));           % Quasi-Steady a
        a(n+1) = a_qs + (a(n)-a_qs)*exp(-dt/tau);   % Update a because of lagging
    end
    
    px_tot = sum(p_x,2);
    py_tot = sum(p_y,2);
    pt_tot = sum(p_t,2);
end