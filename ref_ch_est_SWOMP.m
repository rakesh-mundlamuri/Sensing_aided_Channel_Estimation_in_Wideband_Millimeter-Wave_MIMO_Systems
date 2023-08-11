function [estimated_Hk] = ref_ch_est_SWOMP(r,Nr,Nfft,Psi_at_ar,var_n,txdataF)

epsilon = var_n;
for i = 1:Nfft
    %Gamma(:,:,i) = squeeze(Phi(:,:,i))*Psi_at_ar;  % sw-omp
    Gamma(:,:,i) = txdataF(i)*Psi_at_ar;  % sw-omp
end
estimated_h = zeros(size(Psi_at_ar,2),Nfft);
p_swomp = [];
MSE =0;
%%%%%%%%%%sw-omp%%%%%%%%%%%%%%%%%%%%%%%%%%
y_actual = r;
y_residual = y_actual;

while(1)
    for k = 1:Nfft
        c(:,k)= squeeze(Gamma(:,:,k))'*y_residual(:,k);
    end
    c = abs(c);
    [~,Index] = max(sum(c,2));
    p_swomp = [p_swomp Index];
    
    x_upt_temp = [];
    for k = 1:Nfft
        x_upt_temp(:,k) = pinv(squeeze(Gamma(:, p_swomp,k))) * y_actual(:,k);
        y_residual(:,k) = y_actual(:,k) - squeeze(Gamma(:,p_swomp,k))*x_upt_temp(:,k);
        MSE = MSE + y_residual(:,k)'*y_residual(:,k);
    end
    MSE = MSE/(Nfft*Nr);
    if(MSE<epsilon)
        break;
        
    end
end

for ind = 1 : length(p_swomp)
    estimated_h(p_swomp(ind),:)=x_upt_temp(ind,:);
end

estimated_Hk = Psi_at_ar*estimated_h;
end

