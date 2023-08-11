function [ch_est] = radar_sbl_algo(y,omega,Dictsize_delay,L_prime,Psi,Nfft)

alpha_init = pinv(omega)*y;
gamma = ones(L_prime*Dictsize_delay,1);
Gamma = diag(gamma);
inv_Gamma = inv(Gamma);
z = omega*alpha_init;
zeta = size(y,1)/norm(y-z)^2;
alpha_prev = alpha_init;
k=0;
mse = zeros(1000,1);

while(1)
    k=k+1;
    % output node update
    inv_Gamma_omega_H = inv_Gamma*omega';
    sigma_y_t = (1/zeta)*eye(length(y)) + omega*inv_Gamma_omega_H;
    sigma_t = inv_Gamma - ((inv_Gamma_omega_H/sigma_y_t)*omega*inv_Gamma);
    alpha_est = zeta*(sigma_t)*omega'*y;

    % hyper parameters update
    gamma = real(1./(abs(alpha_est).^2 + diag(sigma_t)));
    inv_Gamma = inv(diag(gamma));
    z = omega*alpha_est;    
    zeta = size(y,1)/norm(y-z)^2;
    mse(k) = norm(alpha_est - alpha_prev)^2/norm(alpha_prev)^2;
    if(mse(k)<1e-8)
        break;
    end
    if(k==1000)
        break;
    end
    alpha_prev = alpha_est;
end

ch_est = zeros(size(Psi,1),size(Psi,3));
for i = 1:Nfft
    ch_est(:,i) = squeeze(Psi(:,:,i))*alpha_est;
end

end
