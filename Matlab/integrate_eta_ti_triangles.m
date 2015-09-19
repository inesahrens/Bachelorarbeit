function result = integrate_eta_ti_total(edges, sur, m, n )
%integrate_eta_ti_total berechnet2 den Vektor int eta t_i dx für alle i
%sur sind alle umgebenden dreiecke also eine n+1*m+1 x 6 matrix
    intEtaTi = integrate_eta_ti_triangles(eta, edges, m, n ); 
    
    result = zeros((m+1)*(n+1),1); 
    for i=1:(m+1)*(n+1)
        result(i) = surroundings_not_null(intEtaTi,sur,i,1,3) + surroundings_not_null(intEtaTi,sur,i,2,3) + surroundings_not_null(intEtaTi,sur,i,3,2) + surroundings_not_null(intEtaTi,sur,i,4,2) + surroundings_not_null(intEtaTi,sur,i,5,1) + surroundings_not_null(intEtaTi,sur,i,6,1);  
    end


end

function result = integrate_eta_ti_triangles(eta, edges, m, n )
%integrate_eta_ti_total berechnet für jedes Dreieck int eta t_i für alle Ti
%die an dem Dreieck anliegen. Für die anderen Ti wäre das Integral sowieso
%0. Daraus ergibt sich, dass wir eine 2*n*m x 3 matrix haben
    
    result = zeros(2*n*m,3) ;
    
    
    result(:,1) = 1/24 * (2 * eta(edges(:,1)) + eta(edges(:,2)) + eta(edges(:,3))) ; 
    result(:,2) = 1/24 * (eta(edges(:,1)) + 2 * eta(edges(:,2)) + eta(edges(:,3))) ; 
    result(:,3) = 1/24 * (eta(edges(:,1)) + eta(edges(:,2)) + 2 * eta(edges(:,3))); 
end
