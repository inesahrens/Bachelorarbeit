function result = G_2(v, eta, c, v0, m, n )
%G2 berechnet G2(v,eta)= -c(v-v0) falls -c(v-v0) <= eta
%                        eta      falls -cv < eta <  -c(v-v0) 
%                        -cv      falls  eta <=  -cv 

    result = zeros((n+1)*(m+1),1); 
    
    for i=1:(n+1)*(m+1)
        if ( -c*(v(i)-v0(i)) <= eta(i) ) 
            result(i) = -c*(v(i)-v0(i)) ; 
        elseif ( - c * v(i) < eta(i) && eta(i) <  -c*(v(i)-v0(i)) )
            result(i) = eta(i); 
        elseif ( eta(i) <= -c* v(i) )
            result(i) = - c * v(i) ; 
        end
    end
    
end

