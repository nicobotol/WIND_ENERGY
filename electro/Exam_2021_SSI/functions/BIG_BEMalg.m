function [cp, cT,MB_normal,MB_tangential, p_normal, p_tangential, ...
            T_normal, T_tangential] = ...
               BIG_BEMalg(R,B,rho,V0,lambda,theta,CL,CD,CM,aoa_ser,glc, ...
                    r,c,beta,tc)
    %% Description
    % This function implements the BEM algorithm
    %% General Information
    % Version: 1
    % Date: 20/09/2021
    % Authors: Sowmya, Philipp, Carlos
    % Denmark Technical University (DTU)
    % Wind Turbine Technologies and Aerodynamics
    % Assignment 1
    %% Function dictionary
    % _____________________________________________________________________
    % INPUTS
    % - r       ---> Radial blade element position [m]
    % - R       ---> Radius of the rotor [m]
    % - B       ---> Number of blades [-]
    % - rho     ---> Air density [kg/m3]
    % - V       ---> Mean wind speed [m/s]
    % - omega   ---> Rotational speed [rad/s]
    % - theta   ---> Pitch angle [deg]
    % - beta    ---> Local twist angle [deg] 
    % - c       ---> Chord of the blade element airfoil [m]
    % - CL      ---> Array with lift coefficient [-]
    % - CD      ---> Array with drag coefficient [-] 
    % - CM      ---> Array with moment coefficient [-]
    % - aoa_ser ---> Series of aoa for which CL,CD, and CM are given [deg]
    % - glc     ---> Glauert correction flag [1 · else]
    % _____________________________________________________________________
    % OUTPUTS
    % - a                  ---> Induced factor normal [-]
    % - aprime             ---> Induced factor tangential [-]
    % - pt                 ---> Local tangential force [N/m]
    % - pn                 ---> Local normal force [N/m]
    % - F                  ---> Correction factor [-]
    % _____________________________________________________________________
    % AUXILIARY
    % - A             ---> Auxiliary matrix for storing the data [-]
    % - count         ---> Auxiliary variable to count BEM iterations [-]
    % - eps           ---> Auxiliary variable for convergence [-]
    % - eps_lim       ---> Tolerance limit for BEM convergence [-]
    % - phi           ---> Flowangle [deg]
    
    % _____________________________________________________________________
    % *********************************************************************
    % Algorithm ===========================================================
    % Step 1 - Initialize a and aprime
    % Step 2 - Compute Flow Angle $Phi$ from Equ. 6.7
    %          Extra: Compute Factor F (Prandtl)
    % Step 3 - Compute local AoA from Equ. 6.6
    % Step 4 - Read C_l($alpha$) and C_d($alpha$) from table look-up
    % Step 5 - Compute C_n and C_t from Equ. 6.12 and 6.13
    % Step 6 - Calculate a and aprime from Equ. 6.23 and 6.24
    %          Prandtl's Tip Loss Factor: Use Equ. 6.35 and 6.36
    % Step 7 - IF a and aprime has changed more than tolerance, go to Step 2,
    %          ELSE Finish
    % Step 8 - Integrate over Blades for resulting Thrust and PowerThis Algo
    % =====================================================================
    %% Parameters
    % *********************************************************************
    eps_lim = 10e-10; % Set a tolerance limit for the convergence
    tc_tolerance = 1e-5; % Tolerance to choose the different tc values
    Nitmax = 1000;
    % *********************************************************************
    %% Operation
    % *********************************************************************
    
NE = size(r,1); % Number of blade elements

    PNplot = zeros((NE+1),1);
    PTplot = zeros((NE+1),1);
    a_plot = zeros(NE,1);
    aprime_plot = zeros(NE,1);
    F_plot = zeros(NE,1);
    
    tc_different = get_different(tc, tc_tolerance);
    
    
    mbending_normal = zeros(NE+1,1);
    mbending_tangential = zeros(NE+1,1);


    for ii = 1:NE % Loop over all the blade elements
        
        tc_be   = tc(ii);
        % Get the thickness of the current blade element
        idx_tc  = which_tc(tc_different,tc_be);
        
        Cl(1,:) = CL(:,idx_tc);
        Cd(1,:) = CD(:,idx_tc);
        Cm(1,:) = CM(:,idx_tc);
        
        r_e     = r(ii);    % NEEDS TO BE IN DICTIONARY
        c_e     = c(ii);    % NEEDS TO BE IN DICTIONARY
        beta_e  = beta(ii); % NEEDS TO BE IN DICTIONARY

        % Here it starts the original BEM
        omega = lambda*V0/R;
        
        % Initialize variables
        a       = 0;                   % Step 1 % Induced Factor normal
        aprime  = 0;                   % Induced Factor tangential
        eps     = 1;                   % Error (has to be larger than eps_lim)
        count   = 0;                   % Iteration counter
        
        while eps > eps_lim && count <= Nitmax % While the convergence is not achieved
            
            % Step 2
            % Calculate the flow angle
            phi     = atand(((1-a)*V0) / ((1+aprime)*omega*r_e));
            
            % Extra (Prandtls Tip Loss)
            % empirical value, NOT acosd. Also include abs(phi)
            F       = (2/pi)*acos(exp(-(B/2)*((R-r_e)/(r_e*abs(sind(phi))))));
            
            % Step 3
            % Compute the local angle of attack
            aoa     = phi - beta_e - theta;
            
            % Step 4 given Cl and Cd
            sigma   = (c_e*B)/(2*pi*r_e);
            [Cl_af, Cd_af, ~] = get_lift_and(Cl, Cd, Cm, aoa_ser,aoa);
            
            Cn      = Cl_af*cosd(phi) + Cd_af*sind(phi);                      % Step 5 normal force
            Ct      = Cl_af*sind(phi) - Cd_af*cosd(phi);                      % tangential force
            
            if glc == 0 || abs(count) <= eps_lim
                
                CT  = ((1-a).^2 * Cn*sigma)/((sind(phi)).^2);
                % Step 6
                an      = 1 / (( (4*(sind(phi)).^2) / (sigma*Cn) ) +1);
                %disp(count)
                
            else
                
                if a <= 1/3
                    an = 1/(((4*F*(sind(phi)).^2)/(sigma*Cn))+1);
                    CT = 4*an*(1-an)*F;
                else
                    astar = CT/(4*F*(1-0.25*(5-3*a)*a));
                    an = 0.1*astar + (1-0.1)*a;
                    %CT = 4*an*(1-0.25*(5-3*an)*an)*F;
                    CT = (1-a)^2*Cn*sigma / (sind(phi))^2;
                    
                end
                
            end
            
            aprime2 = 1 / ( ( (4*F*sind(phi)*cosd(phi)) / (sigma*Ct) ) -1);
            
            eps     = max(abs(a - an),abs(aprime - aprime2));           % Step 7
            a       = an;
            aprime  = aprime2;
            
            
            count   = count + 1;
        end
        
        Vrel    = V0*(1-a) / sind(phi);                             % Step 8
        pN      = 0.5 * rho * Vrel^2 * c_e * Cn;
        pT      = 0.5 * rho * Vrel^2 * c_e * Ct;
        
        %% Compute the loads
        
        PNplot(ii) = pN;
        PTplot(ii) = pT;



        a_plot(ii) = a;
        aprime_plot(ii) = aprime;
        F_plot(ii) = F;

        %%% This is for question 4 · Bending moments at the root
        mbending_normal(ii) = PNplot(ii) * (r(ii) - r(1));
        mbending_tangential(ii) = PTplot(ii) * (r(ii) - r(1));
    end
    
    PNplot(NE+1) = 0;
    PTplot(NE+1) = 0;

    p_normal     = PNplot(1:NE);
    p_tangential = PTplot(1:NE);
    
    r_full = [r' R];
    power  = lambda*V0/R*B*trapz(r_full, r_full.*PTplot');
    cp     = power/(0.5*rho*(V0^3)*pi*(R^2));
    thrust = B * trapz(r_full,PNplot);
    cT     = thrust/ (0.5*rho*V0^2*pi*R^2);
    
    %% This is for Q4
    MB_normal     = trapz(r_full,mbending_normal);
    MB_tangential = trapz(r_full,mbending_tangential);

    %% Shear forces at the root
    T_normal      = trapz(r_full, PNplot);
    T_tangential  = trapz(r_full, PTplot);
end