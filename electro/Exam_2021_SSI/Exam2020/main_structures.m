clear; close all; clc;

%% Exercise 3 - Static Deflections

                    % ########################
                    V0 = 11; % VARY FILE HERE! 6, 11 or 20
                    % ########################
                    
%% Inputs
                    
bladedata   = importdata('bladestruc.txt');
x           = bladedata(:,1);                      % radius [m]
v           = bladedata(:,2);                      % structural pitch [deg]
mass        = bladedata(:,3);                      % distributed mass [kg/m]
EI1         = bladedata(:,4);                      % moment of inertia around PI 1 [Nm2]
EI2         = bladedata(:,5);                      % moment of inertia around PI 2 [Nm2]
twist       = bladedata(:,6);                      % twist angle [deg]

structural_data = "bladestruc.txt";
data_struc  = importdata(structural_data);

NE = size(data_struc,1);

N           = size(x,1);

if V0 == 6
    pitch = 0.896;              % deg
    omega = 0.6283;             % rad/s
    loadsfile = 'loads6.txt';
    
elseif V0 == 11
    pitch = 0;                  % deg
    omega = 0.9253;             % rad/s
    loadsfile = 'loads11.txt';
    
elseif V0 == 20
    pitch = 17.35;              % deg
    omega = 1.0053;             % rad/s
    loadsfile = 'loads20.txt';
else
    warning('ERROR!')
end

loads   = importdata(loadsfile); clear loadsfile;

py      = loads.data(:,3); py = py*1000; % N
pz      = loads.data(:,2); pz = pz*1000; % N

%% Pre-Allocation
Ty      = zeros(N,1);
Tz      = zeros(N,1);
My      = zeros(N,1);
Mz      = zeros(N,1);
M1      = zeros(N,1);
M2      = zeros(N,1);
k1      = zeros(N,1);
k2      = zeros(N,1);
ky      = zeros(N,1);
kz      = zeros(N,1);
th_y    = zeros(N,1);
th_z    = zeros(N,1);
uy      = zeros(N,1);
uz      = zeros(N,1);

%% Step 1: Internal Forces

for ii = N:-1:2
    Ty(ii-1) = Ty(ii) + 0.5*(py(ii-1)+py(ii))*(x(ii)-x(ii-1));
    Tz(ii-1) = Tz(ii) + 0.5*(pz(ii-1)+pz(ii))*(x(ii)-x(ii-1));
    
    My(ii-1) = My(ii) - Tz(ii) * (x(ii)-x(ii-1)) - ...
              ((1/6)*pz(ii-1)+(1/3)*pz(ii)) * (x(ii)-x(ii-1))^2;
    Mz(ii-1) = Mz(ii) + Ty(ii) * (x(ii)-x(ii-1)) + ...
              ((1/6)*py(ii-1)+(1/3)*py(ii)) * (x(ii)-x(ii-1))^2;
end


%% Step 2: Curv(i)ature

for ii = 1:N
    M1(ii)   = My(ii) * cosd(twist(ii)+v(ii)+pitch) - Mz(ii) * sind(twist(ii)+v(ii)+pitch);
    M2(ii)   = My(ii) * sind(twist(ii)+v(ii)+pitch) + Mz(ii) * cosd(twist(ii)+v(ii)+pitch);
    
    k1(ii)   = M1(ii) / (EI1(ii));
    k2(ii)   = M2(ii) / (EI2(ii));
    
    kz(ii)   = -k1(ii)*sind(twist(ii)+v(ii)+pitch) + k2(ii)*cosd(twist(ii)+v(ii)+pitch);
    ky(ii)   = k1(ii)*cosd(twist(ii)+v(ii)+pitch) + k2(ii)*sind(twist(ii)+v(ii)+pitch);
end



%% Step 3: Angle and Deflection

for ii = 1:N-1
    th_y(ii+1) = th_y(ii) + 0.5*(ky(ii+1)+ky(ii))*(x(ii+1)-x(ii));
    th_z(ii+1) = th_z(ii) + 0.5*(kz(ii+1)+kz(ii))*(x(ii+1)-x(ii));
    
    uy(ii+1)   = uy(ii) + th_z(ii)*(x(ii+1)-x(ii)) + ...
                ((1/6)*kz(ii+1)+(1/3)*kz(ii)) * (x(ii+1)-x(ii))^2;
    uz(ii+1)   = uz(ii) - th_y(ii)*(x(ii+1)-x(ii)) - ...
                ((1/6)*ky(ii+1)+(1/3)*ky(ii)) * (x(ii+1)-x(ii))^2;
end

k1 = rad2deg(k1);
k2 = rad2deg(k2);
kz = rad2deg(kz);
ky = rad2deg(ky);
%% Plots

plotter(1,x,py,pz,'p_y','p_z','p(x) [N]','x [m]','Northwest')
plotter(2,x,Ty,Tz,'T_y','T_z','T(x) [N]','x [m]','Northeast')
plotter(3,x,My,Mz,'M_y','M_z','M(x) [Nm]','x [m]','Southeast')
plotter(4,x,ky,kz,'k_y','k_z','k(x) [deg/m]','x [m]','Southwest')
plotter(5,x,th_y,th_z,'theta_y','theta_z','theta(x) [deg]','x [m]','Southwest')
plotter(6,x,uy,uz,'u_y','u_z','u(x) [m]','x [m]','Northwest')

%% Outputs
output_name         = 'statics_out.mat';
delete(output_name);

statics_OUT.x       = x;
statics_OUT.py      = py;
statics_OUT.pz      = pz;
statics_OUT.Ty      = Ty;
statics_OUT.Tz      = Tz;
statics_OUT.My      = My;
statics_OUT.Mz      = Mz;
statics_OUT.ky      = ky;
statics_OUT.kz      = kz;
statics_OUT.th_y    = th_y;
statics_OUT.th_z    = th_z;
statics_OUT.uy      = uy;
statics_OUT.uz      = uz;

save(output_name);
        
%% Static deflection and modeshapes

count = 1;
for ii = 1:2:2*(N)
    m(ii,ii) = 1.1*mass(count);
    m(ii+1,ii+1) = 1.1*mass(count);

    count = count+1;
end

F = zeros(2*N,2*N);    % Flexibility Matrix


for ii = 1:2*N
    
    p       = zeros(2*N,1);
    p(ii)   = 1;
    py      = p(1:2:end);
    pz      = p(2:2:end);
    
    statics_OUT = statics_fun(bladedata,py,pz,0);
    
    uy = statics_OUT.uy;
    uz = statics_OUT.uz;
    
    u            = zeros(2*N,1);
    u(1:2:end)   = uy;
    u(2:2:end)   = uz;
    
    F(:,ii)      = u;
end

[V,D] = eigs(F*m);
omega_1 = sqrt(1/D(1,1));
omega_2 = sqrt(1/D(2,2));
omega_3 = sqrt(1/D(3,3));


%% Mode shapes

MS1_y = V(1:2:end,1) /max(abs(V(:,1)));      %  1 st mode shape y direction
MS1_z = V(2:2:end,1) /max(abs(V(:,1)));      %  1 st mode shape z direction
MS2_y = V(1:2:end,2) /max(abs(V(:,2)));      %  2 nd mode shape y direction
MS2_z = V(2:2:end,2) /max(abs(V(:,2)));      %  2 nd mode shape z direction
MS3_y = V(1:2:end,3) /max(abs(V(:,3)));      %  3 rd mode shape y direction
MS3_z = V(2:2:end,3) /max(abs(V(:,3)));      %  3 rd mode shape z direction

figure('Name','Mode Shapes')
subplot(3,3,1:3)
    plot(x,MS1_y,'DisplayName','y mode')
    hold on
    plot(x,MS1_z,'--','DisplayName','z mode')
    title('1^{st} Mode Shape Flapwise')
    legend('show','Location','Northwest')
    xlim([0 100])
    
subplot(3,3,4:6)
    plot(x,MS2_y,'DisplayName','y mode')
    hold on
    plot(x,MS2_z,'--','DisplayName','z mode')
    title('1^{st} Mode Shape Edgewise')
    legend('show','Location','Southwest')
    xlim([0 100])
    
subplot(3,3,7:9)
    plot(x,MS3_y,'DisplayName','y mode')
    hold on
    plot(x,MS3_z,'--','DisplayName','z mode')
    title('2^{nd} Mode Shape Flapwise')
    legend('show','Location','Northwest')
    xlim([0 100])

