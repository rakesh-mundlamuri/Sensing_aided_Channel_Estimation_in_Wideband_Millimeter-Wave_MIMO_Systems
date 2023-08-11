function [estimated_angles] = SWOMP_algo_angle(r,Nr,pil_len,var_n,txdataF,err_angle_search_range,DictSize_angle,AoA)

epsilon = var_n;

for ang_err = 1
    %[Psi_sbl_angle,dict_angles_sbl]        = Psi_gen_sbl_angle(Nr,AoA,err_angle_search_range/ang_err,DictSize_angle);
    [Psi_sbl_angle,dict_angles_sbl]        = Psi_gen_sbl_angle(Nr,AoA,err_angle_search_range,DictSize_angle);
    Psi_at_ar = Psi_sbl_angle;
    dict_angles = dict_angles_sbl;
    Gamma=[];
    for i = 1:pil_len
        Gamma(:,:,i) = txdataF(i)*Psi_at_ar;  % sw-omp
    end
    %estimated_h = zeros(size(Psi_at_ar,2),pil_len);
    p_swomp = [];
    MSE =0;
    %%%%%%%%%%sw-omp%%%%%%%%%%%%%%%%%%%%%%%%%%
    y_actual = r;
    y_residual = y_actual;

    while(1)
        c = [];
        for k = 1:pil_len
            c(:,k)= squeeze(Gamma(:,:,k))'*y_residual(:,k);
        end
        c = abs(c);
        [~,Index] = max(sum(c,2));
        p_swomp = [p_swomp Index];

        x_upt_temp = [];
        for k = 1:pil_len
            x_upt_temp(:,k) = pinv(squeeze(Gamma(:, p_swomp,k))) * y_actual(:,k);
            y_residual(:,k) = y_actual(:,k) - squeeze(Gamma(:,p_swomp,k))*x_upt_temp(:,k);
            MSE = MSE + y_residual(:,k)'*y_residual(:,k);
        end
        MSE = MSE/(pil_len*Nr);
        %     if(length(p_swomp)==num_paths)
        %         break;
        if(MSE<epsilon)
            break;

        end
    end
    estimated_angles = dict_angles(sort(p_swomp));
    AoA = estimated_angles;
end

end

