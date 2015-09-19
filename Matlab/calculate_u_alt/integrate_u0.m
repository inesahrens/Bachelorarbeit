function links = integrate_u0(edges, allSurroundingTriangles, integral_v, u0, m, n )

    links = zeros((m+1)*(n+1),1);
    for i=0:m
        % berechnung des linken Randes
        surrounders = allSurroundingTriangles(i*(n+1)+1,:);
        intvsurround = integral_of_v_surroundings(surrounders, integral_v); 
        links(i*(n+1)+1) = integrate_u0_one_point_left(surrounders, edges, intvsurround, u0 ); 
%         links(i*(n+1)+1) = boundary(i*(n+1)+1, allSurroundingTriangles, integral_v, edges, u0)
        % berechnung des rechten Randes
        surrounders2 = allSurroundingTriangles(i*(n+1)+n+1,:);
        intvsurround2 = integral_of_v_surroundings(surrounders2, integral_v); 
        links(i*(n+1)+n+1) = integrate_u0_one_point_right(surrounders2, edges, intvsurround2, u0 ) ; 
%         links(i*(n+1)+n+1) = boundary(i*(n+1)+n+1, allSurroundingTriangles, integral_v, edges, u0)
        % berechnung des linken zweiten Randes
        surrounders3 = allSurroundingTriangles(i*(n+1)+2,:); 
        intvsurround3 = integral_of_v_surroundings(surrounders3, integral_v); 
        links(i*(n+1)+2) = integrate_u0_one_point_left2(surrounders3, edges, intvsurround3, u0 ) ;
%         links(i*(n+1)+2) = boundary(i*(n+1)+2, allSurroundingTriangles, integral_v, edges, u0)
        % berechnung des rechten zweiten Randes 
        surrounders4 = allSurroundingTriangles(i*(n+1)+n,:);
        intvsurround4= integral_of_v_surroundings(surrounders4, integral_v); 
        links(i*(n+1)+n) = integrate_u0_one_point_right2(surrounders4, edges, intvsurround4, u0 ) ;
%         links(i*(n+1)+n) = boundary(i*(n+1)+n, allSurroundingTriangles, integral_v, edges, u0)
    end
    
    % diese Werte sollen 0 gesetzt werden, damit die Randwerte beachtet
    % werden. 
    for i=0:m
        links(i*(n+1)+1) = 0;
        links(i*(n+1)+1+n) = 0;
     end
end
% 
% function result = boundary(j, allSurroundingTriangles, integral_v, edges, u0)
%     surrounders = allSurroundingTriangles(j,:);
%     intvsurround = integral_of_v_surroundings(surrounders, integral_v); 
%     result = integrate_u0_one_point_left(surrounders, edges, intvsurround, u0 ) ; 
% end

