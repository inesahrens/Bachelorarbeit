function result = integrate_eta_ti_total(edges, sur, eta, m, n )
%INTEGRATE_ETA_TI_TOTAL  berechnet den Vektor \int_\Omega eta T_i dx für alle i
%sur sind alle umgebenden dreiecke also eine n+1*m+1 x 6 matrix

    % berechnet für jedes einzelne Dreieck das Integral
    intEtaTi = integrate_eta_ti_triangles(eta, edges, m, n ); 
    
    result = zeros((m+1)*(n+1),1); 
    % summert alle Integrale der umliegenden Dreiecke des Gitterpunktes i. 
    for i=1:(m+1)*(n+1)
        result(i) = surroundings_not_null(intEtaTi,sur,i,1,3) + surroundings_not_null(intEtaTi,sur,i,2,3) + surroundings_not_null(intEtaTi,sur,i,3,2) + surroundings_not_null(intEtaTi,sur,i,4,2) + surroundings_not_null(intEtaTi,sur,i,5,1) + surroundings_not_null(intEtaTi,sur,i,6,1);  
    end
end

function result = integrate_eta_ti_triangles(eta, edges, m, n )
%INTERATE_ETA_TI_TRIANGLES berechnet für jedes Dreieck \int_\E eta t_i dx
%f+r T_0,T_1 und T_2
% Daraus ergibt sich, dass wir eine 2*n*m x 3 Matrix haben
    
    result = zeros(2*n*m,3) ;
    
    result(:,1) = 1/24 * (2 * eta(edges(:,1)) + eta(edges(:,2)) + eta(edges(:,3))) ; 
    result(:,2) = 1/24 * (eta(edges(:,1)) + 2 * eta(edges(:,2)) + eta(edges(:,3))) ; 
    result(:,3) = 1/24 * (eta(edges(:,1)) + eta(edges(:,2)) + 2 * eta(edges(:,3))); 
end

function result = surroundings_not_null(intEtaTi, sur, i, j, k)
%SURROUNDINGS_NOT_NULL ist eine funktion, die den wert von intEtaTi
%ausgibt, falls sur(i,j) nicht 0 ist. 
        if (sur(i,j)==0)
            result = 0; 
        else
            result=intEtaTi(sur(i,j),k); 
        end

end
