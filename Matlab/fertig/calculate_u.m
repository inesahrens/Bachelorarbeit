function result = calculate_u(edges, allSurroundingTriangles, u0, epsilon, v, m, n)
% CALCULATE_U berechnet u
%   Ziel ist es  
%   - \int_Omega (v^2 + eps) \nabla u0 \nabla T_j 
%   = \sum_i^{(n+1)(m+1)}  \int_Omega (v^2+eps) \nabla u \nabla T_i \nabla T_j
%   für alle j =1,... (m+1)(n+1) zu berechnen
%   Dafür muss die linke und rechte Seite numerisch berechnet werden. 

    % Da für die Berechnung des Integrals für alle Dreiecke das Integral von
    % v^2+eps gebraucht wird, berechnen wir es hier. 
    integral_v = integral_of_v_total(v, edges, epsilon,n,m); 

%   rechte Seite
%   die rechte Seite ist 0 für alle Werte am Rand und sonst die Werte des
%   Integrals. 
    % rechts muss an den rändern 0 gesetzt werden. Dadurch werden die 0
    % Randwerte eingebracht. Da die Matrix sonst singulär wird, muss der
    % Diagonaleintrag 1 gesetzt werden. 
    rechts = integrate_u(allSurroundingTriangles,integral_v, m, n); 
         for i=0:m
            rechts(i*(n+1)+1,: ) = 0;
            rechts(i*(n+1)+1+n,:) = 0;    
            rechts(i*(n+1)+1,i*(n+1)+1) = 1;
            rechts(i*(n+1)+1+n,i*(n+1)+1+n ) = 1;
         end
         
%   linke Seite: - \int_Omega (v^2 + eps) \nabla u0 \nabla T_j = - \int_Omega (v^2 + eps) \nabla T_i \nabla T_j * u_0 

    links = - rechts * u0; 
    % links muss an den Rändern 0 gesetzt werden. Dadurch werden die 0 Randwerte eingebracht.  
    for i=0:m
       links(i*(n+1)+1) = 0;
       links(i*(n+1)+1+n) = 0;
    end         
    
    % hier wird x aus  Ax=b berechnet. 

    result = rechts\links; 

    result = result + u0; 

end

