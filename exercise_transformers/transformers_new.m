clear all
close 
clc

%%
freq = 60;
V1_LL = 10e3;
n = 10/0.69;

R1 = 3.8e-3; %
X1 = 22.6e-3;
Z1 = R1 + 1j*X1;

R2 = 0.042e-3;
X2 = 0.235e-3;
Z2 = R2 + 1j*X2;
Z2_prime = Z2*n^2;

Rm = 132; 
Xm = 432;
Zm = (1/Rm + 1/(1j*Xm))^(-1);

RL = 0.16;
XL = 2*pi*freq*400e-6;
ZL = RL + 1j*XL;
ZL_prime = ZL*n^2;

%% SECONDARY VOLTAGE
V2_LL = V1_LL/n;
disp(strcat('The LL voltage at the secondary is ', num2str(V2_LL), '  [V]'));

%% POWER AT THE PRIMARY
V1_ph = V1_LL/sqrt(3);

ZTOT = Z1 + (1/Zm + 1/(Z2_prime + ZL_prime))^(-1);
I1 = V1_ph / ZTOT;
S1 = 3*V1_ph*conj(I1);
P1 = real(S1);
Q1 = imag(S1);

%% POWER FACTOR OF THE LOAD
VAB = V1_ph - Z1*I1;
VL = VAB*ZL_prime/(Z2_prime + ZL_prime);
IL = VAB / (ZL_prime + Z2_prime);
SL = 3*VL*conj(IL);
PF = real(SL)/abs(SL);

%% EFFICIENCY
eta = real(SL)/real(S1);

%% 
syms ZL2 RL2

ZL2 = RL2*(1 + 1j*0.8);
Ztot = Z1 + (1/Zm + 1/(Z2_prime + ZL2))^(-1);
I1 = V1_ph/Ztot;
Vab = V1_ph - I1*(1/Zm + 1/(Z2_prime + ZL2))^(-1);
Vl = Vab*ZL2/(Z2_prime + ZL2);
Il = Vab/(Z2_prime + ZL2);
Sl = 3*Vl*conj(Il);
eta2 = real(Sl)/real(S1);
eta2 = simplify(eta2);
dif_eta2 = simplify(diff(eta2, RL2));
vpasolve(dif_eta2 == 0, RL2, 200)

