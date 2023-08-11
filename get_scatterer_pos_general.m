function [scatterer_pos,scatterer_pos_cluster] = get_scatterer_pos_general(BS_pos,UE_pos,scatterer_info)

num        = scatterer_info.num_scatterers;
is_cluster = scatterer_info.is_cluster;
num_paths_per_cluster = scatterer_info.num_paths_per_cluster;

max_x = max(BS_pos(1),UE_pos(1));
min_x = min(BS_pos(1),UE_pos(1));

max_y = 10;
min_y = -10;

scatterer_pos(:,1) = min_x + ((max_x-min_x)*rand(num,1));
scatterer_pos(:,2) = min_y + ((max_y-min_y)*rand(num,1));

scatterer_pos_cluster=[];

if (is_cluster==1)
    error_dist = 0.05;
    for i = 1:num
        clus_err = [[0 0];(-error_dist + 2*error_dist*rand(num_paths_per_cluster-1,2))];
        scatterer_pos_cluster = [scatterer_pos_cluster;scatterer_pos(i,:) + clus_err];
    end
    %scatterer_pos = scatterer_pos_clutter;
end

end

