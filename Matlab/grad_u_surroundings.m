function result = grad_u_surroundings(surroundings, gradUTotal)
%grad_u_surroundings berechnet |nable u|^2 von allen umliegenden Dreiecken 
    result = zeros(1,6); 
    for i=1:6
        if (surroundings(i)~=0)
            result(i)=gradUTotal(surroundings(i));
        end
    end
end

