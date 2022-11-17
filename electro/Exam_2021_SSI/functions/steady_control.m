function outControl = steady_control(V_control, V_step, file_pitch, ...
                        pitch_range, pitch_ld, pitch_rd, cp_opt,...
                        lambda_opt, pitch_opt, R, Omega_min, Omega_max, ...
                        rho, P_R, B, CL,CD,CM,aoa_ser,glc, file_bem_data)
    
    pitch_txt = importdata(file_pitch);
    [r, c, beta, tc] = input_bem_data(file_bem_data);
    N_V = size(V_control,2);

    % Calculate the pitch and rotor speed for each wind speed
    
    pitch_lower = zeros(N_V, 1);
    pitch_upper = zeros(N_V, 1);
    
    omega_control = zeros(N_V,1);
    lambda_control = zeros(N_V,1);
    cP_upper = zeros(N_V,1);
    cT_upper = zeros(N_V,1);
    cP_lower = zeros(N_V,1);
    cT_lower = zeros(N_V,1);
    
    Power  = zeros(N_V, 1);
    Thrust = zeros(N_V, 1);
    Power_lo  = zeros(N_V,1);
    Thrust_lo = zeros(N_V,1);
    
    
    for ii = 1:N_V % Loop for the wind speeds
    
        V0          = V_control(ii);

        if V_step == 1
            V0_idx      = find(pitch_txt(:,1) == V0);
            pitch       = pitch_txt(V0_idx,2);
        else
            pitch = interp1(pitch_txt(:,1),pitch_txt(:,2),V0);
        end

        pitch_step  = (pitch_rd + pitch_ld) / (pitch_range-1);
        pitch_bot   = pitch - pitch_ld;
        pitch_top   = pitch + pitch_rd;
        pitch_vals  = pitch_bot:pitch_step:pitch_top;
        N_Pitch     = length(pitch_vals);
    
        % Calculate rotor speed to keep optimum power production
    
        omega_r = lambda_opt * V0 / R;
    
        if omega_r < Omega_min % Mode 1 of operation · Constant rotor speed at
            % the minimum, cp below optimum, pitch kept
            % constant
    
            omega_r = Omega_min;
            lambda0 = omega_r * R / V0;
            pitch0 = pitch_opt;
            pitch_upper(ii) = pitch0;
            pitch_lower(ii) = pitch0;

            [cP,cT,~,~,~,~] = ...
                BIG_BEMalg(R,B,rho,V0,lambda0,pitch0,...
                CL,CD,CM,aoa_ser,glc, r, c, beta, tc);
            cP_up = cP;
            cP_lo = cP;
            cT_up = cT;
            cT_lo = cT;
    
            cP_upper(ii) = cP_up;
            cT_upper(ii) = cT_up;
            cP_lower(ii) = cP_lo;
            cT_lower(ii) = cT_lo;
        elseif omega_r < Omega_max % Mode 2 of operation · Variable rotor speed
            % cp optimum track, pitch kept constant
    
            cP = cp_opt;
            lambda0 = lambda_opt;
            pitch0 = pitch_opt;

            pitch_upper(ii) = pitch0;
            pitch_lower(ii) = pitch0;

            cP_up = cP;
            cP_lo = cP;
            cT_up = cT;
            cT_lo = cT;
    
            cP_upper(ii) = cP_up;
            cT_upper(ii) = cT_up;
            cP_lower(ii) = cP_lo;
            cT_lower(ii) = cT_lo;
    
        else
    
            omega_r = Omega_max;
            lambda0 = omega_r * R / V0;
            pitch0 = pitch_opt;             % We should change this if we want
            % thrust-peak-saving
    
    
            % Calculate the maximum cP to keep the power at rated
            cP_target = P_R / (0.5 * rho * V0^3 * pi * R^2);
    
            [cP,cT,~,~,~,~] = ...
                BIG_BEMalg(R,B,rho,V0,lambda0,pitch0,...
                CL,CD,CM,aoa_ser,glc,  r, c, beta, tc);
    
            if cP > cP_target % We are in mode 4 · Need pitch control to keep
                % power constant
    
                cP_vals = zeros(N_Pitch, 1);
                cT_vals = zeros(N_Pitch, 1);
    
                for jj = 1:N_Pitch
    
                    pitch_0 = pitch_vals(jj);
    
                    [cP_aux, cT_aux,~,~,~,~] = ...
                        BIG_BEMalg(R,B,rho,V0,lambda0,pitch_0,...
                        CL,CD,CM,aoa_ser,glc,  r, c, beta, tc);
                    cP_vals(jj) = cP_aux;
                    cT_vals(jj) = cT_aux;
                end
    
                %           [cP_vals, cT_vals] = clean_cp_ct(cP_vals, cT_vals);
                %           figure('Name', 'CP curves for pitch control');
                %           plot(pitch_vals, cP_vals);
    
                pitch_control   = calculate_pitch(cP_vals, pitch_vals, cP_target, P_R, rho, R, V0);
                pitch_upper(ii) = pitch_control(2);
                pitch_lower(ii) = pitch_control(1);
    
                %%% Get the actual cp and ct for the estimated pitch values
                [cP_up, cT_up,~,~,~,~] = ...
                    BIG_BEMalg(R,B,rho,V0,lambda0,pitch_control(2),...
                    CL,CD,CM,aoa_ser,glc,  r, c, beta, tc);
    
                [cP_lo, cT_lo,~,~,~,~] = ...
                    BIG_BEMalg(R,B,rho,V0,lambda0,pitch_control(1),...
                    CL,CD,CM,aoa_ser,glc,  r, c, beta, tc);
    
    
    
                %%% We actually want pitch towards feather
    
                cP = cP_up;
                cT = cT_up;
    
                cP_upper(ii) = cP_up;
                cT_upper(ii) = cT_up;
                cP_lower(ii) = cP_lo;
                cT_lower(ii) = cT_lo;
    
            else
                cP_up = cP;
                cP_lo = cP;
                cT_up = cT;
                cT_lo = cT;
    
                cP_upper(ii) = cP_up;
                cT_upper(ii) = cT_up;
                cP_lower(ii) = cP_lo;
                cT_lower(ii) = cT_lo;
    
            end
    
    
        end
    
        %%% Calculate the power for the current wind speed
        Power(ii) = 0.5 * rho * pi * R^2 * V0^3 * cP;
        Thrust(ii) = 0.5 * rho * pi * R^2 * V0^2 * cT;
    
        %%% Calculate the power and thrust for the lower pitch case
        Power_lo(ii) = 0.5 * rho * pi * R^2 * V0^3 * cP_lo;
        Thrust_lo(ii) = 0.5 * rho * pi * R^2 * V0^2 * cT_lo;
    
        omega_control(ii) = omega_r;
        lambda_control(ii) = lambda0;
    end
    
    outControl.omega     = omega_control;
    outControl.lambda    = lambda_control;
    outControl.Pow       = Power;
    outControl.Powlo     = Power_lo;
    outControl.Thrust    = Thrust;
    outControl.Thrust_lo = Thrust_lo;
    outControl.cp_lo     = cP_lower;
    outControl.cp_hi     = cP_upper;
    outControl.cT_lo     = cT_lower;
    outControl.cT_hi     = cT_upper;
    outControl.pitch_lo  = pitch_lower;
    outControl.pitch_hi  = pitch_upper;

end