function result = integrate_u0_one_point_right2(surroundingTriangles, edges, intv, u0 )
    % alledges ist eine 6x3 Matrix
    alledges = edges_of_surrounding_triangles(edges, surroundingTriangles); 
    
    result = -1*u0_not_null(alledges(3,3), u0)*intv(3) - u0_not_null(alledges(6,2), u0)*intv(6) ;
end
