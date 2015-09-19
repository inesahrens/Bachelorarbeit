function result = integrate_u0_one_point_left2(surroundingTriangles, edges, intv, u0 )
    % alledges ist eine 6x3 Matrix
    alledges = edges_of_surrounding_triangles(edges, surroundingTriangles)
    a = u0_not_null(alledges(1,2),u0)
    result = - u0_not_null(alledges(1,2),u0)*intv(1) - u0_not_null(alledges(4,1),u0)*intv(4)
end
