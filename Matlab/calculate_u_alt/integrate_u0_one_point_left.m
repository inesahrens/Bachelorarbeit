function result = integrate_u0_one_point_left(surroundingTriangles, edges, intv, u0 )
    % alledges ist eine 6x3 Matrix
    alledges = edges_of_surrounding_triangles(edges, surroundingTriangles); 
    
    result = (2*u0_not_null(alledges(3,2),u0)-u0_not_null(alledges(3,1),u0))*intv(3) + (u0_not_null(alledges(5,1),u0)-u0_not_null(alledges(5,2),u0))*intv(5) + u0_not_null(alledges(6,1),u0)*intv(6);   
end

