function result = integrate_ti(sur, m, n)
%integrate_ti berechnet das Integral ti über Omega
    result = zeros((n+1)*(m+1),1); 
    
    for i=1:(m+1)*(n+1)
       for j=1:6
          if (sur(i,j)~=0)
             sur(i,j)=1/6; 
          end
       end
        result(i) = sum(sur(i,:)); 
    end

end

