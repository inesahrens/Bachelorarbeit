function result = integral_of_v_surroundings(surroundings,integral_v )
    % gibt alle v^2+eps für anliegende Dreiecke aus
    result = zeros(1,6); 
    for i=1:6
        if (surroundings(i)~=0)
            result(i)=integral_v(surroundings(i));
        end
    end
end

