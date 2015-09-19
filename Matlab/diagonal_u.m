function result = diagonal_u(allSurroundingTriangles, integral_v, n, m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    result = zeros((n+1)*(m+1),1) ; 
    
    for i=1:(m+1)*(n+1)
        surrounders = surrounding_triangles_of_a_point(i,allSurroundingTriangles); 
        intvsurround = integral_of_v_surroundings(surrounders, integral_v)
        result(i) =  intvsurround(1) + intvsurround(2) + 2 * intvsurround(3) + 2 * intvsurround(4) + intvsurround(5) + intvsurround(6) ;    
    end
    
end

