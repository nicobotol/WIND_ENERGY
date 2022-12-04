frequency = 60; % frequency [Hz]
omega = 2*pi*frequency; % [rad/s]
n = 10/0.69; % transformer ratio
V1 = 10000; % primary voltage [V]

R1 = 0.042; % primary resistane (seen form the primary) [ohm]
X1 = 0.032; % primary leakage reactance (seen from the primary) [ohm]
Z1 = R1 + 1j*X1; % primary impedence (seen form primary) [ohm]

R2 = 0.04E-3; % secondary resistance (seen form the secondary) [ohm]
X2 = 0.23E-3; % secondary leakage reactance (seen from the secondary) [ohm]
Z2 = R2 + 1j*X2; % impedence of the secondary (seen from the secondary) [ohm]

Xm = 258; % magnetization resistance (seen from primary side) [ohm]
Rm = 786; % magnetization reactance (seen from primary side) [ohm]
Zm = 1/(1/Rm + 1/(1j*Xm)); % magnetization impedence (seen from the primary) [ohm]

Rload = 0.16; % load resistance [ohm]
Lload = 460E-6; % load inductance [H]
Xload = Lload*omega; % load reactance [ohm]

