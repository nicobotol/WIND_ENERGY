%% Exam WTTA 2021
% DTU
% Philipp Janoff
clear; close all; clc;

%% Q1 =====================================================================

thetavals = [0, 90, 180, 270]

for ii = 1:length(thetavals)
    theta       = thetavals(ii);
    V0          = 12;
    R           = 2;
    omega       = 115 * 2*pi/60;
    Wx          = 4;
    Wy          = 0;
    
    % Coordinate Shifting
    x           = -R*sind(theta);
    y           = R*cosd(theta);
    
    % Velocities
    Vrel_x      = omega*y + V0 - Wx;
    Vrel_y      = -omega*x + Wy;
    Vnorm       = (V0-Wx)*sind(theta) - Wy*cosd(theta);
    Vtan        = (V0-Wx)*cosd(theta) + Wy*sind(theta) + omega*R;
    
    Vrel(ii)    = sqrt(Vnorm^2 + Vtan^2);
    alpha(ii)   = rad2deg(atan(Vnorm/Vtan));
end
Vrel
alpha

%% Q2 =====================================================================
T       = 800000;
V0      = 10;
rho     = 1.225;
R       = 70;
A       = pi*R^2;

CT      = T / (0.5*rho*V0^2*A)

syms a
equ2    = 4*a*(1-a) == CT;
sol_a   = double(solve(equ2,a));

a       = real(sol_a(1))
P       = 2*rho*V0^3*a*(1-a)^2*A

%% Q3 =====================================================================

% Cheat in BEM code around. Solves for \approx 10

%% Q4 =====================================================================
R = 30;
theta = 1;
V0 = 9;
u = 6; % u=(1-a)*V0
TSR = 7;
r = 15;
aoa = 2;
beta = 7;

phi = aoa+beta+theta;
omega = TSR*V0/R;
aprime = u/(tand(phi)*omega*r)-1

%% Q5 =====================================================================
F = load('flexibility.mat','F');
F = F.F;    
p = zeros(36,1);
p(1:2:end) = 0;
p(2:2:end) = 1;

for ii = 1:1:5000000
    utest = F*(ii*p);
    utend(ii) = utest(end);
end

const = interp1(utend,1:1:5000000,5)

%% Q6 =====================================================================

% Use main_structures.m and add a factor of 1.1 to mass distribution.


%% Q7 =====================================================================

% derivation done

%% Q8 =====================================================================

Power = [0.7, 2, 3];
Vbot  = [5, 10, 13];
Vtop  = [10, 13, 25];

A = 8;
k = 2;

for ii = 1:length(Power)
    f = exp(-(Vbot(ii)/A)^k)-exp(-(Vtop(ii)/A)^k);
    AEP_ii(ii) = 8760*f*Power(ii)
end

AEP = sum(AEP_ii)

%% Q9 =====================================================================

% General Inputs
f       = 60;               % Hz
w       = 2*pi*f; 
V1      = 10000;            % Volt
atr     = 10000/690;            % ratio

% Resistances
R1      = 3.8e-3; 
X1      = 22.6e-3;
R2p     = atr^2 * 0.042e-3;
X2p     = atr^2 * 0.235e-3;
R3      = 432;
X3      = 132;
RIp     = atr^2 * 0.15;
XIp     = atr^2 * 300e-6 * w; %¿¿¿¿ w ???

% Impedances
Z1 = R1+1j*X1;
Z2_1 = R2p+1j*X2p;

Z2_2 = RIp+1j*XIp;

Z2 = Z2_1 + Z2_2;
Z3 = (R3*1j*X3)/(R3+1j*X3);

Zeq = Z1 + (Z2*Z3)/(Z2+Z3);

%Computations
phi = rad2deg(angle(Z3));

I1  = V1/sqrt(3) / Zeq;     % Use line to line : 1/sqrt(3)
I3  = (V1/sqrt(3)-Z1*I1)/Z3;
I2p = I1-I3;

V2p = V1 - I1*Z1 - Z2_1*I2p;
V2  = 1/atr * V2p;                    % Part a ============================
V2_mag = abs(V2)

S1  = sqrt(3) * V1 * conj(I1);
P1  = real(S1)
Q1  = imag(S1)                      % Part b ==============================

I2  = atr * I2p;
S2  = sqrt(3)*V2*conj(I2);

pf  = real(S2)/abs(S2)              % Part c ==============================

% improve: Add capacitor to load, 
%          thus lower phase between I2 and V2 
%           -> increase real(), decrease imag()

eta = real(S2)/real(S1)             % Part d ==============================
% Power_out / Power_in

