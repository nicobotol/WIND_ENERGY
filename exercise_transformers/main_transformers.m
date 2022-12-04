% Three phase transformer (connected in star) operates at nominal 60Hz and
% has a transformer ratio of line to line voltages 10kV:0.69kV. At the 
% nominal frequency and voltage, the primary resistance and leakage 
% reactance (as measured from primary side) are 4.2 and 32 mOhms. The 
% secondary resistance and leakage reactance (as measured from secondary 
% side) are 0.04 and 0.23 mOhms. The magnetizing reactance and the core 
% resistance (as measured from primary side) are 258 and 786 Ohms.
% The voltage on primary (HV side) is constant 10kV (line to line). The
% three-phase industrial load is connected to secondary (LV side) and is 
% characterized with resistance 0.16 Ohms and inductance 460 uH. Using "T"
% equivalent circuit of the transformer answer following

% load parameters from file
clear 
close all
clc
parameters

%% Voltage at the secondary

V2 = V1 / n;
disp(strcat('The LL voltage at the secondary is ', num2str(V2), '  [V]'))

%% Active and reactive power at HV side

Zload = Rload + 1j*Xload; % impedence of the load [ohm]
Z2_prime = Z2*n^2; % impedence of the secondary repoorted to the primary [ohm]
Zload_prime = Zload*n^2; % impedence of the load seen from the primary [ohm]

Ztot = Z1 + (1/Zm + 1/(Z2_prime + Zload_prime))^(-1);
I1 = V1/sqrt(3)/Ztot; % per-phase current [A]
S1_phase = V1/sqrt(3)*conj(I1); % complex power of the primary (per phase) [VA]
S1_tot = 3*S1_phase; % total complex power of the primary [VA]
P1_tot = real(S1_tot); % total active power of the primary [W]
Q1_tot = imag(S1_tot); % total reactive power of the primary [VAR]
disp(strcat('Active power at the primary', num2str(P1_tot), '[W]'))
disp(strcat('Reative power at the primary', num2str(Q1_tot), '[VAR]'))

%% Power factor of the load
V1_ph = V1/sqrt(3);
Z3 = Z2_prime + Zload_prime;
V2_ph = V1_ph *(1/Z2_prime + 1/Z3)^(-1)/(Z1 + (1/Z2_prime + 1/Z3)^(-1));
I2 = V2_ph / Z3;
VL = V2*Zload_prime/(Z2_prime + Zload_prime);
SL = VL*conj(I2);
PF_load = real(SL)/abs(SL);

%phi_load = atan(X2/R2); % phase angle of the load
%PF_load = cos(phi_load); % power factor of the load
disp(strcat('The power factor of the load is ', num2str(PF_load)))

%% Effiency of the transformer
Zpar = (1/Zm + 1/(Z2_prime + Zload_prime))^(-1);

Vab = V1/sqrt(3)*Zpar/(Zpar + Z1);
I2 = Vab / (Z2_prime + Zload_prime);
S2 = I2*Zload_prime*conj(I2);
P2 = real(S2);
P2_tot = 3*P2;
eta = P2_tot/P1_tot; % efficiency
disp(strcat('The efficiency of the transformer is:', num2str(eta)))

