function [H_time,H_freq,Psi] = gen_channel(Nt,Nr,Nfft,Ncp,Ts,rays,alpha,flag,flag_psi_gen)

rolloff = 0.9;
zt = (0:Nt-1)';
zr = (0:Nr-1)';
num_paths = sum(flag);
num_comm_paths = sum(flag);

% Steering vectors
At = zeros(Nt,num_paths);
Ar = zeros(Nr,num_paths);
for i=1:num_paths
    At(:,i) = exp(1i*pi*cosd(rays.AoD(i))*zt);
    Ar(:,i) = exp(1i*pi*cosd(rays.AoA(i))*zr);
end

H_time = zeros(Nr,Nt,Ncp);
for d=0:Ncp-1
    for paths = 1:num_paths
        H_time(:,:,d+1) = H_time(:,:,d+1) + alpha(paths) ...
            * raised_cosine_filter(d,rays.delay(paths),Ts,rolloff) ...
            * Ar(:,paths) * At(:,paths)';
    end
end %paths

H_freq = zeros(Nr,Nt,Nfft);

for nt = 1:Nt
    for nr = 1:Nr
        H_freq(nr,nt,:) = fft(H_time(nr,nt,:),Nfft);
    end %nr
end %nt

% Psi matrix generation
Psi = [];
if(flag_psi_gen==1)
    beta_k_l = zeros(Nfft,num_paths);
    for k = 0:Nfft-1
        for paths = 1:num_paths
            temp_d = [];
            for d = 0:Ncp-1
                temp_d = [temp_d raised_cosine_filter(d,rays.delay(paths),Ts,rolloff)*exp(-1j*2*pi*d*k/Nfft)];
            end
            beta_k_l(k+1,paths) = sum(temp_d);
        end
    end

%     H_freq = zeros(Nr,Nt,Nfft);
% 
%     for k=0:Nfft-1
%         for paths = 1:num_paths
%             H_freq(:,:,k+1) = H_freq(:,:,k+1) + (alpha(paths)*beta_k_l(k+1,paths)*Ar(:,paths)* At(:,paths)');
%         end
%     end %paths

    ATbar_circ_Ar = zeros(Nt*Nr,num_comm_paths);
    beta_k_l_ext = zeros(Nfft,num_comm_paths);
    ind = 1;
    for i = 1:num_paths
        if(flag(i)==1)
            ATbar_circ_Ar(:,ind) = kron(conj(At(:,i)),Ar(:,i));
            beta_k_l_ext(:,ind)  = beta_k_l(:,i);
            ind=ind+1;
        end
    end

    Psi = zeros(Nt*Nr,num_comm_paths,Nfft);
    for i = 1:Nfft
        Psi(:,:,i) = ATbar_circ_Ar*diag(beta_k_l_ext(i,:));
    end
end

end

