clc;
clear;
close all;

%% Initial Parameters
c                 = physconst('lightspeed');% speed of light in m/s
fc                = 28e9;              % centre frequency
lambda            = c/fc;              % wavelength
SCS               = 120e3;             % subcarrier spacing
Nt_UE             = 1;                 % Number of transmit antennas at the UE
Nr_gNB            = 32;                % Number of receive antennas at the gNB
BS_pos            = [2 3];             % Base station location in xy coordinates
UE_pos            = [7 4];             % UE location in xy coordinates
LOS               = 1;                 % if there is a line of sight set this to 1 or else 0
t_symb            = 1/SCS;
Nfft              = 256;                 % fft size
Ncp               = 34;                  % length of the cyclic prefix of the OFDM symbol in a first slot
sampling_rate     = Nfft*SCS;            % sampling rate
Ts                = 1/sampling_rate;     % sampling period
Gt                = 1;                   % dictionary size for AoD
Gr                = 500;                 % dictionary size for AoA
pil_len_nb        = 16;
pilot_comb        = Nfft/pil_len_nb;
extra_pilots_nb   = 0;
err_angle_range   = 3;                   % angle error in degrees
err_delay_range   = Ts/2;
err_angle_search_range = 2*err_angle_range;
err_delay_search_range = 2*err_delay_range;

SNR_measure       = 0;
DictSize_angle    = 500;
Dictsize_delay    = 50;
NUM_ITER          = 100;
SNR_vec           = -10:5:20;
%SNR_vec           = 20;

% Flags to run different cases
is_cluster            = 0;
num_scatterers        = 9;                 % number of scatterers
num_com_scatterers    = 5;

if (is_cluster == 1)
    str_model               = "cluster_model";
    num_scatterers          = 2;                 % number of scatterers
    num_paths_per_cluster   = 5;
    num_com_scatterers      = num_scatterers;
    num_paths               = num_com_scatterers*num_paths_per_cluster+LOS;
else
    str_model               = "point_scatterers_model";
    num_paths_per_cluster   = 0;
    num_paths               = num_com_scatterers+LOS;
end
scatterer_info.is_cluster            = is_cluster;
scatterer_info.num_scatterers        = num_scatterers;
scatterer_info.num_paths_per_cluster = num_paths_per_cluster;
scatterer_info.num_comm_scatterers   = num_com_scatterers;
scatterer_info.num_paths             = num_paths;

pilot_indices                        = 1:pilot_comb:Nfft;

%% SIMULATION LOOP:
% Psi matrix generation for SWOMP
At_g           = exp(pi*1j * (0:Nt_UE-1)' * (1-2/Gt*(0:Gt-1)));
Ar_g           = exp(pi*1j * (0:Nr_gNB-1)' * (1-2/Gr*(0:Gr-1)));
Psi_at_ar      = kron(conj(At_g),Ar_g);
omega_ls       = pinv(Psi_at_ar'*Psi_at_ar)*Psi_at_ar';

% generate pilots
txdataF = sqrt(1/2)*sqrt(1/Nt_UE)*(sign(randn(Nt_UE,Nfft))+1j*sign(randn(Nt_UE,Nfft)));

% Phi generation
ext_pilots         = txdataF(1,pilot_indices);

% Initialize error variables
err_ch_est_swomp_SNR   = zeros(length(SNR_vec),1);
err_ch_est_SNR_nb      = zeros(length(SNR_vec),1);
err_ch_est_SNR_sbl     = zeros(length(SNR_vec),1);
err_ch_est_SNR_ls      = zeros(length(SNR_vec),1);
SNR_est_avg            = zeros(length(SNR_vec),1);

% Simulation loop over SNR
parfor SNR_dB = 1:length(SNR_vec)
    SNR_vec(SNR_dB)
    err_ch_nb             = zeros(NUM_ITER,1);
    err_ch_sbl            = zeros(NUM_ITER,1);
    err_ch_swomp          = zeros(NUM_ITER,1);
    err_ch_ls             = zeros(NUM_ITER,1);
    SNR_est_collect       = zeros(NUM_ITER,1);

    % Simulation loop over iterations
    for iter = 1:NUM_ITER
        iter
        var_n             = 10^(-SNR_vec(SNR_dB)/10);
        SNR_aux           = zeros(Nfft,1);

        [scatterer_pos,scatterer_pos_cluster]     = get_scatterer_pos_general(BS_pos,UE_pos,scatterer_info);

        if(LOS==1)
            env_pos         = [inf inf;scatterer_pos];
            env_pos_cluster = [inf inf;scatterer_pos_cluster];
            is_LOS          = [1 zeros(1,size(scatterer_pos,1))];
        else
            env_pos         = scatterer_pos;
            env_pos_cluster = scatterer_pos_cluster;
            is_LOS          = zeros(1,size(scatterer_pos,1));
        end

        %% Generate the rays using the created tool
        rays                      = gen_paths(BS_pos,UE_pos,env_pos,is_LOS);
        rays_with_error       = rays;
        rays_with_error.AoA   = rays.AoA + normrnd(0,err_angle_range,1,length(rays.AoA));
        rays_with_error.AoA(rays_with_error.AoA < 0 ) = 180+ rays_with_error.AoA(rays_with_error.AoA < 0 );
        rays_with_error.delay = rays.delay + normrnd(0,err_delay_range,1,length(rays.delay));
        rays_with_error.delay(rays_with_error.delay < 0 ) = 0;

        % generate channel gains
        sigma2            = sort(randfixedsum(num_paths,1,1,0,1),'descend');  % noise variances of each path
        alpha             = zeros(num_paths,1);

        for i = 1:num_paths
            if(rays.is_LOS(i) == 1)
                alpha(i) = normrnd(0,sqrt(sigma2(i)),1)*exp(-1j*2*pi*fc*rays.delay(i));
            else
                alpha(i) = normrnd(0,sqrt(sigma2(i)),1)*exp(-1j*((2*pi*fc*rays.delay(i)) + (2*pi*rand(1,1))));
            end
        end

        [~,H_freq,Psi]                         = gen_channel(Nt_UE,Nr_gNB,Nfft,Ncp,Ts,rays,alpha,ones(num_paths,1),1);

        % Received signal generation
        signal_k          = zeros(Nr_gNB,Nfft);
        for i = 1:Nfft
            signal_k(:,i) = squeeze(H_freq(:,:,i))*txdataF(:,i);
        end

        % Noise generation
        Noise             = sqrt(var_n/2)*(randn(Nr_gNB,Nfft)+1i*randn(Nr_gNB,Nfft));

        % received signal
        r                 = signal_k+Noise;
        ext_r             = r(:,pilot_indices);
        ext_Psi           = Psi(:,:,pilot_indices);

        % Received SNR calculation
        if(SNR_measure == 1)
            for i = 1: Nfft
                pow_sig   = signal_k(:,i)'*signal_k(:,i);
                pow_nn    = Noise(:,i)'*Noise(:,i);
                SNR_aux(i)= pow_sig/pow_nn;
            end
            SNR_est_collect(iter) = 10*log10(mean(SNR_aux));
        end

        ch_est_ls = zeros(Nr_gNB,Nfft);
        for i = 1:Nfft
            alpha_est_ls = txdataF(i)'*omega_ls*r(:,i);
            ch_est_ls(:,i) = Psi_at_ar*alpha_est_ls;
        end

        % channel estimation method: Ideal Sensing Info + LS
        omega_nb       = [];
        for i = 1:pil_len_nb
            omega_nb   = [omega_nb;ext_pilots(1,i)*squeeze(ext_Psi(:,:,i))];
        end
        alpha_est_nb   = pinv(omega_nb)*ext_r(:);
        ch_est_nb      = zeros(Nr_gNB,Nfft);
        for i = 1:Nfft
            ch_est_nb(:,i)     = squeeze(Psi(:,:,i))*alpha_est_nb;
        end

        % channel estimation method: SWOMP-SBL + Sensing Info Error
        est_angles_swomp       = SWOMP_algo_angle(ext_r,Nr_gNB,pil_len_nb,var_n,ext_pilots,err_angle_search_range,DictSize_angle,rays_with_error.AoA);
        L_prime                = length(est_angles_swomp);
        path_ass_est           = zeros(L_prime,1);
        for ang_ind = 1:L_prime
            [~,path_ass_est(ang_ind)] = min(abs(est_angles_swomp(ang_ind)-rays_with_error.AoA));
        end
        [Psi_sbl_delay,dict_delay,angle_append,delay_append,omega_sbl_delay] = Psi_gen_swomp_delay(Nr_gNB,Nfft,Ncp,Ts,rays_with_error,est_angles_swomp,err_delay_search_range,Dictsize_delay,path_ass_est,ext_pilots,pilot_indices);
        ch_est_sbl             = radar_sbl_algo(ext_r(:),omega_sbl_delay,Dictsize_delay,L_prime,Psi_sbl_delay,Nfft);

        % channel estimation method: WB + SWOMP
        ch_est_swomp           = ref_ch_est_SWOMP(r,Nr_gNB,Nfft,Psi_at_ar,var_n,txdataF);

        % NMSE of the channel
        err_ch_ls(iter)        = sum(abs(ch_est_ls(:)-H_freq(:)).^2)/(norm(H_freq(:),'fro')^2);
        err_ch_swomp(iter)     = sum(abs(ch_est_swomp(:)-H_freq(:)).^2)/(norm(H_freq(:),'fro')^2);
        err_ch_nb(iter)        = sum(abs(ch_est_nb(:)-H_freq(:)).^2)/(norm(H_freq(:),'fro')^2);
        err_ch_sbl(iter)       = sum(abs(ch_est_sbl(:)-H_freq(:)).^2)/(norm(H_freq(:),'fro')^2);
        %fprintf('SNR %d dB, iter %d, err_ch_ls %f,err_ch_swomp %f,err_ch_nb %f,err_ch_sbl %f\n',SNR_vec(SNR_dB),iter,err_ch_ls(iter),err_ch_swomp(iter),err_ch_nb(iter),err_ch_sbl(iter));

    end

    err_ch_est_SNR_ls(SNR_dB)          = mean(err_ch_ls);
    err_ch_est_swomp_SNR(SNR_dB)       = mean(err_ch_swomp);
    err_ch_est_SNR_nb(SNR_dB)          = mean(err_ch_nb);
    err_ch_est_SNR_sbl(SNR_dB)         = mean(err_ch_sbl);

    if(SNR_measure==1)
        SNR_est_avg(SNR_dB)            = mean(SNR_est_collect);
    end
end
%% Plots
figure();

plot(SNR_vec,10*log10(err_ch_est_SNR_nb),'-o');
hold on;
grid on;
plot(SNR_vec,10*log10(err_ch_est_swomp_SNR),'->');
plot(SNR_vec,10*log10(err_ch_est_SNR_sbl),'-+');
plot(SNR_vec,10*log10(err_ch_est_SNR_ls),'-x');

title('SNR VS NMSE of channel');
xlabel('SNR(dB)');
ylabel('NMSE(dB)');

str_leg = [];
str_leg = [str_leg ;"NB + Radar info"];
str_leg = [str_leg ;"WB + SWOMP"];
str_leg = [str_leg ;"NB + SBL + Radar info error"];
str_leg = [str_leg ;"Classical LS"];

legend(str_leg);
str_save_fig       = "pilot_comb_"+pilot_comb+"_num_pilots_"+pil_len_nb+"_angle_error_"+err_angle_range+"_ITER_"+NUM_ITER+str_model+".png";
str_save_variables = "pilot_comb_"+pilot_comb+"_num_pilots_"+pil_len_nb+"_angle_error_"+err_angle_range+"_ITER_"+NUM_ITER+str_model+".mat";
saveas(gcf,str_save_fig);

save(str_save_variables,'SNR_vec','err_ch_est_SNR_nb','err_ch_est_swomp_SNR','err_ch_est_SNR_sbl','err_ch_est_SNR_ls');
