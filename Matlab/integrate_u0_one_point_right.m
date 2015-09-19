function result = integrate_u0_one_point_right(surroundingTriangles, edges, intv, u0 )
    % alledges ist eine 6x3 Matrix
    alledges = edges_of_surrounding_triangles(edges, surroundingTriangles); 
    
    result = u0_not_null(alledges(1,3), u0)*intv(1) + (u0_not_null(alledges(2,3),u0)-u0_not_null(alledges(2,2),u0))*intv(2) +  (2*u0_not_null(alledges(4,2),u0)-u0_not_null(alledges(4,3),u0))*intv(4);
end

