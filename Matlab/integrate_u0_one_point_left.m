function result = integrate_u0_one_point_left(alledges, intv, u0 )
    
    result = (2*u0_not_null(alledges(3,2),u0)-u0_not_null(alledges(3,1),u0))*intv(3) + (u0_not_null(alledges(5,1),u0)-u0_not_null(alledges(5,2),u0))*intv(5) + u0_not_null(alledges(6,1),u0)*intv(6);   
end

