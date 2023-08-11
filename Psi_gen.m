function [Psi] = Psi_gen(Nt,Nr,Nfft,Ncp,Ts,rays,flag)
rolloff = 0.9;
zt = (0:Nt-1)';
zr = (0:Nr-1)'; 
num_paths = length(rays.is_LOS);
num_comm_paths = sum(flag);

% Steering vectors        
ATbar_circ_Ar = zeros(Nt*Nr,num_comm_paths);
At = zeros(Nt,num_paths);
Ar = zeros(Nr,num_paths);
ind = 1;
for i=1:num_paths
    if(flag(i)==1)
        At(:,i) = exp(1i*pi*cosd(rays.AoD(i))*zt);
        Ar(:,i) = exp(1i*pi*cosd(rays.AoA(i))*zr);
        ATbar_circ_Ar(:,ind) = kron(conj(At(:,i)),Ar(:,i));
        ind=ind+1;
    end
end

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

% Psi matrix generation
beta_k_l_ext = zeros(Nfft,num_comm_paths);
ind=1;
for i = 1:num_paths
    if(flag(i)==1)        
        beta_k_l_ext(:,ind)  = beta_k_l(:,i);
        ind=ind+1;
    end
end

Psi = zeros(Nt*Nr,num_comm_paths,Nfft);
for i = 1:Nfft
    Psi(:,:,i) = ATbar_circ_Ar*diag(beta_k_l_ext(i,:));
end

end

