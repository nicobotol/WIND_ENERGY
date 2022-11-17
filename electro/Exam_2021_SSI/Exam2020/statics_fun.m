function statics_OUT = statics_fun(bladedata, py, pz, pitch)

% Inputs

x           = bladedata(:,1);                      % radius [m]
v           = bladedata(:,2);                      % structural pitch [deg]
EI1         = bladedata(:,4);                      % moment of inertia around PI 1 [Nm2]
EI2         = bladedata(:,5);                      % moment of inertia around PI 2 [Nm2]
twist       = bladedata(:,6);                      % twist angle [deg]

N           = size(x,1);    % Number of blade points



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

% k1 = rad2deg(k1);
% k2 = rad2deg(k2);
kz = rad2deg(kz);
ky = rad2deg(ky);

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

end