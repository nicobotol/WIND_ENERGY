function theta = find_theta(cp, pitch_vals, cp_max)
    
    AA = cp - cp_max;
    [~,j_min] = min(abs(AA)); %find theta
    
    if j_min+1 > size(AA,1)
        theta = ...
            interp1([cp(j_min-1), cp(j_min)], ...
            [pitch_vals(j_min-1),pitch_vals(j_min)],cp_max);
        
        % If the target point is located between the minimum and
        % the next point
    elseif AA(j_min)*AA(j_min+1) <= 0
        theta = ...
            interp1([cp(j_min), cp(j_min+1)], ...
            [pitch_vals(j_min),pitch_vals(j_min+1)],cp_max);
        
        
        % If the target point is located between the minimum and
        % the previous point
    else
        theta = ...
            interp1([cp(j_min-1), cp(j_min)], ...
            [pitch_vals(j_min-1),pitch_vals(j_min)],cp_max);
        
    end
end