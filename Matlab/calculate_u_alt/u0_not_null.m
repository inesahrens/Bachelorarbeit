function result = u0_not_null(a,u0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
        if (a==0)
            result = 0; 
        else
            result=u0(a); 
        end
end

