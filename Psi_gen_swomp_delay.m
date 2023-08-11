function [Psi,dict_delay,angle_append,delay_append,omega] = Psi_gen_swomp_delay(Nr,Nfft,Ncp,Ts,rays,AoA_known,delay_err,DictSize_delay,path_ass_est_angle_swomp,txdataF,pilot_indices)

rolloff = 0.9;
zr = (0:Nr-1)';
delay_known = rays.delay(path_ass_est_angle_swomp);
dict_delay = [];

for iter = 1:length(path_ass_est_angle_swomp)
    delay_min   = delay_known(iter)-delay_err;
    if(delay_min < 0)
        delay_min = 0;
    end

    delay_max   = delay_known(iter)+delay_err;
    
    delay_step  = (delay_max-delay_min)/(DictSize_delay-1);
    delay_range = delay_min:delay_step:delay_max;
    dict_delay = [dict_delay delay_range];
end

% Steering vectors        
Ar = exp(1i*pi*zr*cosd(AoA_known));

num_paths_delay_dict = length(dict_delay);

beta_d_l = zeros(Ncp,num_paths_delay_dict);
for paths = 1:num_paths_delay_dict
    beta_d_l(:,paths) = raised_cosine_filter(0:Ncp-1,dict_delay(paths),Ts,rolloff);
end
beta_k_l = fft(beta_d_l,Nfft,1);
%beta_k_l_append = repmat(beta_k_l,1,length(AoA_known));
beta_k_l_append = beta_k_l ;
%delay_append = repmat(dict_delay,1,length(AoA_known));
delay_append = 0;

Ar_append = [];
angle_append= [];
% for i =1:length(AoA_known)
%     Ar_append = [Ar_append repmat(Ar(:,i),1,size(beta_k_l,2))];
%     angle_append = [angle_append repmat(AoA_known(i),1,size(beta_k_l,2))];
% end
for i =1:length(AoA_known)
    Ar_append = [Ar_append repmat(Ar(:,i),1,DictSize_delay)];
    angle_append = [angle_append repmat(AoA_known(i),1,DictSize_delay)];
end

% Psi = zeros(Nr,length(AoA_known)*size(beta_k_l,2),Nfft);
% for i = 1:Nfft
%     Psi(:,:,i) = Ar_append*diag(beta_k_l_append(i,:));
% end

Psi = zeros(Nr,length(AoA_known)*DictSize_delay,Nfft);
for i = 1:Nfft
    Psi(:,:,i) = Ar_append*diag(beta_k_l_append(i,:));
end

omega = [];
for i = 1:length(txdataF)
    omega = [omega;txdataF(i)*squeeze(Psi(:,:,pilot_indices(i)))];
end

end
