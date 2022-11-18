line_width = 2;
font_size = 30;

%% Data of the problem
rho = 1.225; % air density (kg/m^3)
Rs = 0.064; % wtg internal resistnace (ohm)
Ls = 0.0018; % wtg internal inductance (H)
flux = 19.49; % rototr magnetic flux (Wb)
poles = 320; % number of poles of the rotor
pp = poles/2; % pole pairs of the rotor
omega_M_max = 30*1.01/pi; % maximum rotational speed (rpm)
R = 89.17; % generator rotor radius (m)
A = pi*R^2; % blade swapt area
cp_opt = 0.465; % optimal cp
lambda_opt = 7.857; % optimal lambda

% transformers parameters
L = 65; % km of cables
c_cable = 0.1E-6; % capacity for unit length (F)
l_cable = 0.5E-3; % inductance for unit length (H)
r_cable = 0.1; % resistance for unit length (ohm)
Cc = c_cable*L; % cable capacity
Lc = l_cable*L;
Rc = r_cable*L;
n = 4/33; % transformer ratio
Lm = 0.046; % transformer magnetizing inductance (H)
L1 = 0.6E-3; % inductance of the primary (H)
L2_prime = L1; % inductance of the secondary referred to the primary
R1 = 0.025; % resitance of the primary (Ohm)
R2_prime = 0.025; % resitance of the secondary referred to the primary (Ohm)
f_grid = 50; % grid side frequency (Hz)
Vb_ll = 4000; % line to line voltage VSC-B (V)
Vpoc = 33000; % line to line voltage POC (V)