function links = integrate_u0(edges, allSurroundingTriangles, integral_v, u0, m, n )

    links = zeros((m+1)*(n+1),1);
    for i=0:m
        % berechnung des linken Randes
        [intvsurround, alledges] = surrounding_things(allSurroundingTriangles, integral_v, edges, i*(n+1)+1); 
        links(i*(n+1)+1) = integrate_u0_one_point_left(alledges, intvsurround, u0 ) ; 
        % berechnung des rechten Randes
        [intvsurround2, alledges2] = surrounding_things(allSurroundingTriangles, integral_v, edges, i*(n+1)+1+n); 
        links((i+1)*(n+1)) = integrate_u0_one_point_right(alledges2, intvsurround2, u0 ) ; 
        % berechnung des linken zweiten Randes
        [intvsurround3, alledges3] = surrounding_things(allSurroundingTriangles, integral_v, edges, i*(n+1)+2); 
        links(i*(n+1)+2) = integrate_u0_one_point_left2(alledges3, intvsurround3, u0 ); 
        % berechnung des rechten zweiten Randes  
        [intvsurround4, alledges4] = surrounding_things(allSurroundingTriangles, integral_v, edges, i*(n+1)+n); 
        links((i+1)*(n+1)-1) = integrate_u0_one_point_right2(alledges4, intvsurround4, u0 ); 
    end
    
    % diese Werte sollen 0 gesetzt werden, damit die Randwerte beachtet
    % werden. 
%     for i=0:m
%         links(i*(n+1)+1) = 0;
%         links(i*(n+1)+1+n) = 0;
%      end
end

function [intvsurround, alledges] = surrounding_things(allSurroundingTriangles, integral_v, edges, i)
        surrounders = allSurroundingTriangles(i,:);
        intvsurround = integral_of_v_surroundings(surrounders, integral_v);
        alledges = edges_of_surrounding_triangles(edges, surrounders);
end


function alledges = edges_of_surrounding_triangles(edges, surroundings)
% gibt zu einem Punkt i alle Eckpunkte der umliegenden Dreiecke an. 
    alledges = zeros(6,3);
    for j=1:6
       if (surroundings(j) ~= 0)
            alledges(j,:)= edges(surroundings(j),:);
       end
    end

end

function result = integrate_u0_one_point_left(alledges, intv, u0 )    
    result = (2*u0_not_null(alledges(3,2),u0)-u0_not_null(alledges(3,1),u0))*intv(3) + (u0_not_null(alledges(5,1),u0)-u0_not_null(alledges(5,2),u0))*intv(5) + u0_not_null(alledges(6,1),u0)*intv(6);   
end

function result = integrate_u0_one_point_right(alledges, intv, u0 )    
    result = u0_not_null(alledges(1,3), u0)*intv(1) + (u0_not_null(alledges(2,3),u0)-u0_not_null(alledges(2,2),u0))*intv(2) +  (2*u0_not_null(alledges(4,2),u0)-u0_not_null(alledges(4,3),u0))*intv(4);
end

function result = integrate_u0_one_point_left2(alledges, intv, u0 )
    result = -1* u0_not_null(alledges(1,2),u0)*intv(1) - u0_not_null(alledges(4,1),u0)*intv(4);   
end

function result = integrate_u0_one_point_right2(alledges, intv, u0 )    
    result = -1*u0_not_null(alledges(3,3), u0)*intv(3) - u0_not_null(alledges(6,2), u0)*intv(6) ;
end

function result = u0_not_null(a,u0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
        if (a==0)
            result = 0; 
        else
            result=u0(a); 
        end
end



