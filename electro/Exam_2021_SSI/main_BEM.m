%%% _______________________________________________________________________
%%% BEM Main Script
%%% Assignment 1 ---- Wind Turbine Technologies and Aerodynamics 46300
%%% Authors:
%%%  - Sowmya
%%%  - Philipp
%%%  - Carlos

clear; close all; clc;

% _________________________________________________________________________
% · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · 
%% General Parameters
% Conversion factors ******************************************************
rpm2rads = 2 * pi / 60;      % Conversion factor rpm to rad/s [-]
% Atmospheric conditions **************************************************
rho      = 1.225;            % Atmospheric density [kg/m3]
% Wind turbine definition *************************************************
B        = 3;                % Number of blades [-]
P_R      = 10.5e6;             % Rated power of the wind turbine [W]
R        = 89.17;            % Rotor radius [m]
TS_max   = 90;               % Maximum tip speed allowed [m/s]
                             % (typically 90 for offshore, 70 for onshore)
Omega_min = 6.0 * rpm2rads;  % Minimum rotor speed [rad/s]
Omega_max = TS_max / R;      % Maximum rotor speed [rad/s]
V_ci     = 4;                % Cut-in Wind speed [m/s]
V_co     = 25;               % Cut-out Wind speed [m/s]
V_rated  = 11.4;             % Rated wind speed of the wind turbine [m/s]
% BEM resolution **********************************************************
glc      = 1;                % Glauert correction flag for the BEM code

% _________________________________________________________________________
% · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % ·
%% Importing files and data

functions_folder = ".\functions";
input_folder     = ".\inputs";

addpath(functions_folder);
addpath(input_folder);

file_aoa_ser    = "aoa_series.dat";
aoa_ser    = importdata(file_aoa_ser);

file_bem_data   = "bladedat.txt";
Data = readtable(file_bem_data);
NE = size(Data,1);
clearvars Data;

file_blade_aero = "aerodynamics_bem.dat";
positions = [1 2 3];

file_performance = "BEM_performance.mat";

% The aerodynamic files should be in increasing distance to hub
files_aero = ["cylinder.txt" "FFA-W3-600.txt" "FFA-W3-480.txt" ...
    "FFA-W3-360.txt" "FFA-W3-301.txt" "FFA-W3-241.txt"];
tc_files     = [100, 60, 48, 36, 30.1, 24.1]; % Thicknesses ratios

table_aerodynamics = importdata(file_blade_aero);
[CL,CD,CM] = extract_from_matrix(table_aerodynamics, positions);

[r, c, beta, tc] = input_bem_data(file_bem_data);

% _________________________________________________________________________
% · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · 
%% Parameters Question 1 · Find optimum cP

calculate_cp_opt = 0;       % Flag to indicate wether is required to 
                            % calculate the optimum cp

% Values of angle of attack for the airfoils interpolation
interpdata.alpha_min    = -180;    % Minimum value of angle of attack 
interpdata.alpha_max    = 180;     % Maximum value of angle of attack
interpdata.d_alpha      = 0.1;     % Distance between aoa

% Values for calculating the curves
% Values of the blade pitch (according to design data)
curvesdat.pitch_min    = -2;
curvesdat.pitch_max    = 5;
curvesdat.d_pitch      = 0.1;
% Values of the TSR (according to design data)
curvesdat.lambda_min   = 5;
curvesdat.lambda_max   = 10;
curvesdat.d_lambda     = 0.1;


% _________________________________________________________________________
% · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · 
%% Parameters Question 2 · Running at optimum cP

% Keep it at 1 until the functions are updated for variable step size
V_step  = 1;

V_control = V_ci:V_step:V_co;
N_V = size(V_control,2);

% _________________________________________________________________________
% · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · 
%% Parameters Question 3 · Pitch control

file_pitch = 'Pitch.txt';   % File with orientative pitch angles
% Pitch Inputs to create the extended Array based on TABLE pitch value
pitch_range = 41;           % Number of pitch points
pitch_ld    = 40;           % Extension to the left side (for 
pitch_rd    = 2;

% _________________________________________________________________________
% · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · 
%% Parameters Exercise 3 · Structural calculations
file_structural  = 'bladestruc.txt';
file_out_statics = 'output_statics.mat';

% _________________________________________________________________________
% · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · 
%% Parameters Question 4 · Bending moments
% Compute the bending moments with respect to the root of the blade root
% for a series of wind speeds
%VQ4 = [4 8 11 18 25]; % Wind speeds to calculate the bending moments
VQ4 = [6 11 20];


% _________________________________________________________________________
% · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · % · 
%% Parameters Question 5
% Calculate average energy production and difference with early cut-out
A = 9;         % Weibull Scaling factor of the site 
k = 1.9;       % Weibull Shape factor of the site 
V_co_new = 20; % New cut-out wind speed to compare with the actual one

%% Import data

load BEM_performance.mat;           % Load the MAT file with the 

% pitch_txt = importdata(file_pitch);

cp_opt      = BEM_out.cpopt;       % Extract the optimum power coefficient
lambda_opt  = BEM_out.tsropt;      % Optimum tip-speed-ratio
pitch_opt   = BEM_out.pitchopt;    % Optimum pitch angle
cp_curves   = BEM_out.cp;          % Cp curves for question 1
cT_curves   = BEM_out.cT;          % CT curves for question 1
lam_vals    = BEM_out.lambdavals;  % Lambda values for the curves
the_vals    = BEM_out.thetavals;   % Pitch values for the curves


V_omegaN   = TS_max / lambda_opt;
V_omegamin = Omega_min * R / lambda_opt;

outControl = steady_control(V_control, V_step, file_pitch, ...
                        pitch_range, pitch_ld, pitch_rd, cp_opt,...
                        lambda_opt, pitch_opt, R, Omega_min, Omega_max, ...
                        rho, P_R, B, CL,CD,CM,aoa_ser,glc, file_bem_data);

omega_control   = outControl.omega;
lambda_control  = outControl.lambda;
Power           = outControl.Pow;
Power_lo        = outControl.Powlo;
Thrust          = outControl.Thrust;
Thrust_lo       = outControl.Thrust_lo;
cP_lower        = outControl.cp_lo;
cP_upper        = outControl.cp_hi;
cT_lower        = outControl.cT_lo;
cT_upper        = outControl.cT_hi;
pitch_lower     = outControl.pitch_lo;
pitch_upper     = outControl.pitch_hi;

%% PLOTTING 
figure('Name', 'Pitch control')
plot(V_control, pitch_upper, 'k-o');
hold on;
plot(V_control, pitch_lower,'k--');
grid on
legend ('High Pitch Values','Low Pitch Values', 'Interpreter', 'Latex',...
    'Location', 'northwest');
xlabel('$V$ [m/s]', 'Interpreter', 'Latex')
ylabel('Pitch $\theta_{P}$ [\textordmasculine]', 'Interpreter','latex');
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('pitch_control','-dpdf');

figure('Name', 'Power curve');
plot(V_control, Power*1e-6,'k-o');
hold on
plot(V_control,Power_lo*1e-6,'k--');
grid on
legend ('High Pitch Values','Low Pitch Values', 'Interpreter', 'Latex',...
    'Location', 'northwest');
xlabel('$V$ [m/s]', 'Interpreter', 'Latex')
ylabel('Power [MW]', 'Interpreter','latex');
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('power_curve_control','-dpdf');

figure('Name', 'Thrust curve');
plot(V_control, Thrust*1e-6,'k-o');
hold on
plot(V_control,Thrust_lo*1e-6,'k--');
grid on
legend ('High Pitch Values','Low Pitch Values', 'Interpreter', 'Latex',...
    'Location', 'northwest');
xlabel('$V$ [m/s]', 'Interpreter', 'Latex')
ylabel('Thrust [MN]', 'Interpreter','latex');
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('thrust_curve_control','-dpdf');

%% TILES

figure('Name', 'Comparing curves');
compare_curves = tiledlayout(6,1);

ax1 = nexttile;
plot(V_control, omega_control / rpm2rads,'k-', ...
    'DisplayName','Rotor speed')
ylabel('Rotor speed [RPM]', 'Interpreter','latex');
hold on
yline(Omega_min / rpm2rads, 'k--', ...
    'DisplayName','Minimum rotor speed')
yline(Omega_max / rpm2rads, 'k-.', ...
    'DisplayName','Minimum rotor speed')
xlim([V_ci, V_co]);
ylim([ (min(omega_control / rpm2rads) - max(omega_control / rpm2rads*0.1)) ...
       (max(omega_control / rpm2rads)  + max(omega_control / rpm2rads*0.1)) ]);
grid on
xline(V_rated,'k--')
xline(V_omegaN,'k--')
xline(V_omegamin,'k--')

ax2 = nexttile;
plot(V_control, Power/(10^6),'k-', 'DisplayName','WT Power')
xlim([V_ci, V_co]);
ylabel('Power [MW]', 'Interpreter','latex');
ylim([ (min(Power/(10^6)) - max(Power/(10^6)*0.1)) ...
       (max(Power/(10^6))  + max(Power/(10^6)*0.1)) ]);
grid on
hold on
xline(V_rated,'k--')
xline(V_omegaN,'k--')
xline(V_omegamin,'k--')

ax3 = nexttile;
plot(V_control, cP_upper,'k-', 'DisplayName','C_{P}');
xlim([V_ci, V_co]);
ylim([ (min(cP_upper) - max(cP_upper*0.1)) ...
       (max(cP_upper)  + max(cP_upper*0.1)) ]);
ylabel('$C_{P}$', 'Interpreter','latex');
grid on
hold on
xline(V_rated,'k--')
xline(V_omegaN,'k--')
xline(V_omegamin,'k--')


ax4 = nexttile;
plot(V_control, pitch_upper,'k-');
xlim([V_ci, V_co]);
ylim([ (min(pitch_upper) - max(pitch_upper*0.1)) ...
       (max(pitch_upper)  + max(pitch_upper*0.1)) ]);
ylabel('Pitch [deg]', 'Interpreter','latex');
grid on
hold on
xline(V_rated,'k--')
xline(V_omegaN,'k--')
xline(V_omegamin,'k--')


ax5 = nexttile;
plot(V_control, Thrust/(10^3),'k-');
xlim([V_ci, V_co]);
ylim([ (min(Thrust/(10^3)) - max(Thrust/(10^3))*0.1) ...
       (max(Thrust/(10^3))  + max(Thrust/(10^3))*0.1) ]);
ylabel('Thrust [MN]','Interpreter','latex');
grid on
hold on
xline(V_rated,'k--')
xline(V_omegaN,'k--')
xline(V_omegamin,'k--')

ax6 = nexttile;
plot(V_control, cT_upper,'k-');
xlim([V_ci, V_co]);
ylim([ (min(cT_upper) - max(cT_upper*0.1)) ...
       (max(cT_upper)  + max(cT_upper*0.1)) ]);
ylabel('$C_{T}$','Interpreter','latex');
grid on
hold on
xline(V_rated,'k--')
xline(V_omegaN,'k--')
xline(V_omegamin,'k--')

xticklabels(ax1,{});
xticklabels(ax2,{});
xticklabels(ax3,{});
xticklabels(ax4,{});
xticklabels(ax5,{});

xlabel('Wind Speed [m/s]','Interpreter','latex');
% title(compare_curves, ...
%     'Characteristic curves of the DTU-10MW Reference Wind Turbine', ...
%     'Interpreter', 'Latex')

compare_curves.TileSpacing ='compact';
print('tiles','-dpdf','-fillpage');

%%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%%
%% Question 4 · Bending moments
%%% _______________________________________________________________________

NVQ4 = length(VQ4);
Normal_bending     = zeros(NVQ4,1);
Tangential_bending = zeros(NVQ4,1);
Normal_force       = zeros(NE,NVQ4);
Tangential_force   = zeros(NE,NVQ4);

for ii = 1:NVQ4
    
   V0 = VQ4(ii); % Take the ii wind speed
   
   % Find the index of the wind speed in the series of TSR and Pitch, to
   % extract the values
   idx_vel     = find_me_the_index(V_control, V0);
   lambda0     = lambda_control(idx_vel);
   pitch0      = pitch_upper(idx_vel);
   omega0      = omega_control(idx_vel);

   % Call BEM to get the normal bending moment and the tangential bending
   % moment
   [~, ~, MB1, MB2, PF1, PF2] = ...
                BIG_BEMalg(R,B,rho,V0,lambda0,pitch0,...
                CL,CD,CM,aoa_ser,glc,  r, c, beta, tc);
    
    % Store the bending moments in their respective arrays
    Normal_bending(ii)     = MB1;
    Tangential_bending(ii) = MB2;
    % Store the normal and tangential forces in their respective arrays
    Normal_force(:,ii)       = PF1;
    Tangential_force(:,ii)   = PF2;

end


%% Plot the bending moments
figure('Name', 'Bending moments')
plot(VQ4,Normal_bending.*10^(-6),'k-o', ...
    'DisplayName', 'Bending moment flapwise'); 
hold on; 
plot(VQ4,Tangential_bending.*10^(-6),'k-d', ...
    'DisplayName', 'Bending moment edgewise');
xlabel('Wind speed [m/s]', 'Interpreter', 'Latex');
ylabel('Moment [$MN \cdot m$]', 'Interpreter', 'Latex');
grid on
legend('Flapwise Bending Moment','Edgewise Bending Moment');
%flapwise is normal bending moment, edgewise is tangential bending moment;
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('bending_moments','-dpdf');

%%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%% ··· %%%
%% Question 5 · Energy production
%%% _______________________________________________________________________

AEP = average_energy_production(Power,V_control,k,A);
fprintf('The AEP with cut-out at %i m/s is %4.2f GWh\n', V_co, ...
    AEP * 10^(-9));
% Plot the Weibull distribution
v_weibull = 0:0.01:30;
f_v = weibull(v_weibull, k, A);
figure('Name', 'Weibull plot');
plot(v_weibull,f_v, 'k-');
xlabel('Wind speed [m/s]', 'Interpreter', 'Latex');
ylabel('Probability f(v)', 'Interpreter', 'Latex');
title('\textbf{Weibull probability distribution}','Interpreter', 'Latex');
grid on;
text(15,0.09,{['Scaling factor, A = ',num2str(A)],...
    ['Shape factor, k = ', num2str(k)]}, 'Interpreter', 'Latex')
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('weibull_distribution','-dpdf');

%% What if the cut-out wind speed changes?

idx_vel = find_me_the_index(V_control, V_co_new);

% Vectors with the power curve for early co
P_earlyco = Power(1:idx_vel); 
V_earlyco = V_control(1:idx_vel);

% Call function to calculate the AEP in the early cut-out
AEP_early = average_energy_production(P_earlyco,V_earlyco,k,A);

AEP_Difference = AEP - AEP_early;
AEP_Diff_Relative = AEP_Difference / AEP * 100; % Percentage of the losses

fprintf('The AEP with early cut-out at %i m/s is %4.2f GWh\n', V_co_new, ...
    AEP_early * 10^(-9));
fprintf('The difference is %4.2f GWh\n', AEP_Difference * 10^(-9));
fprintf('And the relative loss is %4.2f%%\n', AEP_Diff_Relative);

%% Plots of the comparison of cut-out with Weibull distribution

idx_velci = find_me_the_index(v_weibull, V_ci);
idx_velco = find_me_the_index(v_weibull, V_co);
idx_velco_new = find_me_the_index(v_weibull, V_co_new);
f_used = f_v;
f_used_new = f_v;
f_used(1:idx_velci) = 0;
f_used_new(1:idx_velci) = 0;
f_used(idx_velco:end) = 0;
f_used_new(idx_velco_new:end) = 0;


figure('Name', 'Weibull distribution with early cut-out');
area(v_weibull,f_used,'FaceColor', '#A9A9A9', ...
    'DisplayName','Standar flow','LineStyle',':');
hold on;
area(v_weibull,f_used_new,'FaceColor', '#708090', ...
    'DisplayName','Early cut-out flow','LineStyle',':');
plot(v_weibull,f_v,'k-','DisplayName','Wind distribution');
xlabel('Wind speed [m/s]', 'Interpreter', 'Latex');
ylabel('Probability f(v)', 'Interpreter', 'Latex');
title_string = {'\textbf{Weibull probability distribution}', ...
    ['Used wind flow for early cut-out at $V_{co} = ',...
    num2str(V_co_new), ' m/s$']};
title(title_string, 'interpreter', 'latex', 'fontsize', 12);
grid on;
text(20,0.05,{['Scaling factor, A = ',num2str(A)],...
    ['Shape factor, k = ', num2str(k)]}, 'Interpreter', 'Latex')
legend show
set(gcf, 'PaperPosition', [0 0 15 10]); 
set(gcf, 'PaperSize', [15 10]); 
print('weibull_distribution_filled','-dpdf');

