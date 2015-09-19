function result = G_2_eta(v, eta, c, v0, m, n )
%G2eta berechnet G2eta(v,eta)= 0          falls -c(v-v0) < eta oder eta < -cv 
%                         1          falls -cv < eta <  -c(v-v0) 
%                         in [0,1]   falls -c(v-v0) = eta oder eta = -cv 

    result = zeros((n+1)*(m+1),1); 
    
    for i=1:(n+1)*(m+1)
        if ( -c*(v(i)-v0(i)) < eta(i) || eta(i) < -c* v(i)) 
            result(i) = 0 ; 
        elseif ( - c * v(i) < eta(i) && eta(i) <  -c*(v(i)-v0(i)) )
            result(i) = 1; 
        elseif ( -c*(v(i)-v0(i)) == eta(i) || eta(i) == -c* v(i) )
            %hier immer 0.5. kann sein, dass ich das später ändern muss. 
            result(i) = 0.5 ; 
        end
    end
    
end
