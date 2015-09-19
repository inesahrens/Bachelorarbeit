function alledges = edges_of_surrounding_triangles(edges, surroundings)
% gibt zu einem Punkt i alle Eckpunkte der umliegenden Dreiecke an. 
    alledges = zeros(6,3);
    for j=1:6
       if (surroundings(j) ~= 0)
            alledges(j,:)= edges(surroundings(j),:);
       end
    end

end

