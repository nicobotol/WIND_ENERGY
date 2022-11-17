%% WIND TURBINE TECHNOLOGY AND AERODYNAMICS
%Exercise 4 Vertical Axis Wind Turbine

clear; close all; clc;


%% VAWT INPUTS
dt      = 0.001;     % Timestep
ntime   = 10000;
nblades = [1,2,3];  % Number of blades



%% Initialize and Loop

for ii = 1:length(nblades)
    nb = nblades(ii);
    [t(:,ii),px(:,ii),py(:,ii),pt(:,ii),CP(:,ii),CT(:,ii)] = VAWT(nb,dt,ntime);
end

mean_cp = mean(CP)
mean_ct = mean(CT)


%% Plot
figure('Name','Power & Thrust Coefficients')                 % CP CT PLOT
    subplot(3,3,1:3)                                % total
        plot(t(:,1),CP(:,1),'k-','DisplayName','C_P')
        hold on
        plot(t(:,1),CT(:,1),'b-','DisplayName','C_T')
        xlabel('Time t [s]')
        ylabel('Coefficients')
        legend('show','Location','Southwest')
        title('Coefficients: B=1')
    
    subplot(3,3,4:6)                                % total
        plot(t(:,2),CP(:,2),'k-','DisplayName','C_P')
        hold on
        plot(t(:,2),CT(:,2),'b-','DisplayName','C_T')
        xlabel('Time t [s]')
        ylabel('Coefficients')
        legend('show','Location','Southwest')
        title('Coefficients: B=2')
        
    subplot(3,3,7:9)                                % total
        plot(t(:,3),CP(:,3),'k-','DisplayName','C_P')
        hold on
        plot(t(:,3),CT(:,3),'b-','DisplayName','C_T')
        xlabel('Time t [s]')
        ylabel('Coefficients')
        legend('show','Location','Southwest')
        title('Coefficients: B=3')
            
        
    

                
        
figure('Name','Total Loads')
    subplot(3,3,1:3)
        plot(t(:,1),px(:,1),'b-','DisplayName','p_{x,tot}');
        hold on
        plot(t(:,1),py(:,1),'r-','DisplayName','p_{y,tot}');
        legend('show')
        title('B=1')
        xlabel('Time t [s]')
        ylabel('total loads [Nm]')
        
    subplot(3,3,4:6)
        plot(t(:,2),px(:,2),'b-','DisplayName','p_{x,tot}')
        hold on
        plot(t(:,2),py(:,2),'r-','DisplayName','p_{y,tot}')
        legend('show')
        title('B=2')
        xlabel('Time t [s]')
        ylabel('total loads [Nm]')
        
    subplot(3,3,7:9)
        plot(t(:,3),px(:,3),'b-','DisplayName','p_{x,tot}')
        hold on
        plot(t(:,3),py(:,3),'r-','DisplayName','p_{y,tot}')
        legend('show')
        title('B=3')
        xlabel('Time t [s]')
        ylabel('total loads [Nm]')
        
        
        
        
% azi_idx = theta(:,1);
% azi_idx(azi_idx<2*pi) = 0;
% azi_idx = find(azi_idx,1,'first')-1;
%         
% figure('Name','AoA over Azimuth Position')
%     plot(azi(1:azi_idx,1),rad2deg(alpha(1:azi_idx,1)),'k-','DisplayName','Blade 1')
%     hold on
%     %plot(azi(1:azi_idx,2),rad2deg(alpha(1:azi_idx,2)),'b--','DisplayName','Blade 2')
%     %plot(azi(1:azi_idx,3),rad2deg(alpha(1:azi_idx,3)),'r-.','DisplayName','Blade 3')
%     set(gca,'XTick',0:pi/2:2*pi) 
%     set(gca,'XTickLabel',{'0','\pi/2','\pi','3*\pi/2','2*\pi'})
%     ylabel('Angle of Attack [deg]')
%     xlabel('Azimuth Position [rad]')
    
    
%% Function

function [t,px_tot,py_tot,pt_tot,CP,CT] = VAWT(B,dt,ntime)

    %% Data Input
    file_name = 'airfoil.txt';
    airfoil   = importdata(file_name);

    aoa         = deg2rad(airfoil.data(:,1));           % AoA given in DEG, then to RAD
    lift_coeff  = airfoil.data(:,2);
    drag_coeff  = airfoil.data(:,3);

    % Parameters
    R       = 3;        % Radius                     [m]
    V0      = 8;        % Freestream Velocity        [m/s]
    omega   = 14;       % Angular Velocity           [rad/s]
    S       = 0.2;      % Surface Area               [Bc/R]
%     c       = 0.1;      % Chord                      [m]
    rho     = 1.225;    % Density                    [kg/m3]
    tau     = 2*R/V0;   % Time Relaxation constant   [s]

    a       = zeros(ntime,1);
    a(1)    = 0; % Initialize a=0. Redundant but just for protocol.

    for n = 1:ntime
        
        t(n)        = (n-1)*dt;
        theta1      = omega*t(n);
        W_global    = a(n) * V0;
        
        for ii = 1:B                                         % for each blade
            theta(n,ii) = theta1 + 2*pi*(ii-1)/B; % Shifting Theta per Blade
            azi(n,ii)   = mod(theta(n,ii),2*pi);
            Wx          = W_global * (1 - 0.3*sin(theta(n,ii)));
            Wy          = 0.3*Wx*cos(theta(n,ii));
            
            % Coordinate Shifting
            x           = -R*sin(theta(n,ii));
            y           = R*cos(theta(n,ii));
            
            % Velocities
            Vrel_x      = omega*y + V0 - Wx;
            Vrel_y      = -omega*x + Wy;
            Vnorm       = (V0-Wx)*sin(theta(n,ii)) - Wy*cos(theta(n,ii));
            Vtan        = (V0-Wx)*cos(theta(n,ii)) + Wy*sin(theta(n,ii)) + omega*R;
            Vrel        = sqrt(Vnorm^2 + Vtan^2);
            
            alpha(n,ii) = atan(Vnorm/Vtan);
            pitch       = 0;                                % NO PITCH!
            phi         = alpha(n,ii) + pitch;
            
            % Interpolate
            Cl(n,ii)    = interp1(aoa,lift_coeff,alpha(n,ii));
            Cd(n,ii)    = interp1(aoa,drag_coeff,alpha(n,ii));
            
            c           = S*R/B;
            
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