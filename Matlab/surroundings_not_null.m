function result = surroundings_not_null(intEtaTi, sur, i, j, k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        if (sur(i,j)==0)
            result = 0; 
        else
            result=intEtaTi(sur(i,j),k); 
        end

end

