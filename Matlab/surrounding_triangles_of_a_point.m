function surroundings = surrounding_triangles_of_a_point(i,surroundingTriangles )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
    %Berechnet zu einem Punkt i, alle umliegenden Dreiecke
    surroundings = zeros(1,6);
    surroundings(:) = surroundingTriangles(i,:); 
end

