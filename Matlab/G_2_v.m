function result = G_2_v(v, eta, c, v0, m, n )
%G2v berechnet G2v(v,eta)= -c         falls -c(v-v0) < eta oder eta < -cv 
%                         0          falls -cv < eta <  -c(v-v0) 
%                         in [-c,0]  falls -c(v-v0) = eta oder eta = -cv 

    result = zeros((n+1)*(m+1),1); 
    
    for i=1:(n+1)*(m+1)
        if ( -c*(v(i)-v0(i)) < eta(i) || eta(i) < -c* v(i)) 
            result(i) = -c ; 
        elseif ( - c * v(i) < eta(i) && eta(i) <  -c*(v(i)-v0(i)) )
            result(i) = 0; 
        elseif ( -c*(v(i)-v0(i)) == eta(i) || eta(i) == -c* v(i) )
            %hier immer -c/2 kann sein, dass ich das später ändern muss. 
            result(i) = -c/2 ; 
        end
    end
    
end