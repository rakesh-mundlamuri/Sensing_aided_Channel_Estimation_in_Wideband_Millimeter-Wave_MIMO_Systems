function [Psi,dict_angles] = Psi_gen_sbl_angle(Nr,AoA,angle_err,DictSize_angle)

zr = (0:Nr-1)';
AoA_known = AoA;
dict_angles = [];

for iter = 1:length(AoA)
    if(Nr>1)
        theta_AoA_min   = AoA_known(iter)-angle_err;
        if(theta_AoA_min<0)
            theta_AoA_min = 0;
        end

        theta_AoA_max   = AoA_known(iter)+angle_err;
        if(theta_AoA_max>180)
            theta_AoA_max = 180;
        end

        theta_AoA_step  = (theta_AoA_max-theta_AoA_min)/(DictSize_angle-1);
        theta_AoA_range = theta_AoA_min:theta_AoA_step:theta_AoA_max;
    else
        theta_AoA_range = 0;
    end
    dict_angles = [dict_angles theta_AoA_range];
end

% Steering vectors        
Ar = exp(1i*pi*zr*cosd(dict_angles));

Psi = Ar;

end

