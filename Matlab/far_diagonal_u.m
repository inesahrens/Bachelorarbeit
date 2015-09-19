function result = far_diagonal_u(allSurroundingTriangles,integral_v, n, m)

    result = zeros((n+1)*(m+1),1); 
    
    for i=1:(n+1)*(m+1)
       surrounders = surrounding_triangles_of_a_point(i,allSurroundingTriangles);
       intvsurround = integral_of_v_surroundings(surrounders, integral_v); 
       result(i) = - intvsurround(4) - intvsurround(5);    
    end


end

