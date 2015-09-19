function v = u0_vektor(m,n, v,u0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



for i=1:(m+1)*(n+1) 
    if (mod(i,n+1)==0 || mod(i,n+1)==1)
        v(i)=u0(i); 
    end
end

end

