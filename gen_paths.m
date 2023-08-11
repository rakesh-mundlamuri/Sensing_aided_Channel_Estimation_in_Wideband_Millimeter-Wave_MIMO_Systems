function path = gen_paths(BS_pos,UE_pos,scatterer_pos,is_LOS)
c = 3e8; % speed of light in m/s
% angles are measured in anticlockwise direction wrt to x axis

for i = 1:length(is_LOS)
    if(is_LOS(i))
        AoD = atand((BS_pos(2)-UE_pos(2))/(BS_pos(1)-UE_pos(1)));
        AoA = atand((UE_pos(2)-BS_pos(2))/(UE_pos(1)-BS_pos(1)));
        distance = sqrt((BS_pos(1)-UE_pos(1))^2+(BS_pos(2)-UE_pos(2))^2);
        delay = distance/c;
        path.delay(i) = delay;
        path.distance(i) = distance;
        path.is_LOS(i) = is_LOS(i);
        if(AoA<0)
            path.AoA(i) = AoA+180;
        else
            path.AoA(i) = AoA;
        end
        if(AoD<0)
            path.AoD(i) = AoD+180;
        else
            path.AoD(i) = AoD;
        end
    end
    if(is_LOS(i) == 0)
        AoD = atand((BS_pos(2)-scatterer_pos(i,2))/(BS_pos(1)-scatterer_pos(i,1)));
        AoA = atand((UE_pos(2)-scatterer_pos(i,2))/(UE_pos(1)-scatterer_pos(i,1)));
        distance = sqrt((BS_pos(1)-scatterer_pos(i,1))^2+(BS_pos(2)-scatterer_pos(i,2))^2)+sqrt((UE_pos(1)-scatterer_pos(i,1))^2+(UE_pos(2)-scatterer_pos(i,2))^2);
        delay = distance/c;
        path.delay(i) = delay;
        path.distance(i) = distance;
        path.is_LOS(i) = is_LOS(i);
        if(AoA<0)
            path.AoA(i) = AoA+180;
        else
            path.AoA(i) = AoA;
        end
        if(AoD<0)
            path.AoD(i) = AoD+180;
        else
            path.AoD(i) = AoD;
        end
    end
end

end

