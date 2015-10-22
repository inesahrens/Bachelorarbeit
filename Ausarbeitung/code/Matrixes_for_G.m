function [A, B, D] =  Matrixes_for_G(u, edges, allSurroundingTriangles,...
    m, n )
%Berechnet alle Matrizen und Vektoren, die für die Berechnung von G
%gebraucht werden, außer e. 

    % grad_u_total berechnet |\nabla u|^2 für alle Dreiecke der 
    % Triangulierung. 
    gradUTotal = grad_u_total(u, edges); 
    
    % A = \int\Omega |\nabla u|^2 T_i T_j dx
    A = integrate_u_ti_tj(gradUTotal, allSurroundingTriangles, m, n); 

    % B = \int\Omega \nabla T_i \nabla T_j dx
    B = integrate_gti_gtj(allSurroundingTriangles, m, n); 

    % D = \int\Omega T_i T_j dx
    D = integrate_ti_tj(allSurroundingTriangles, m, n); 
end

function gradU = grad_u_total(u, edges)
%GRAD_U_TOTAL berechnet |\nabla u|^2 für alle Dreiecke der Triangulierung. 
% Dabei ist  |\nabla u|^2 = (u31-u21)^2 + (u31-u21)^2 + (u32-u22)^2 +
% (u32-u22)^2 Die Formel kann in der Bachelorarbeit nachgeschlagen werden. 
    gradU = ( u(edges(:,3),1) - u(edges(:,2),1) ).^2 ...
        + ( u(edges(:,1),1) - u(edges(:,2),1) ).^2 ...
        + ( u(edges(:,3),2) - u(edges(:,2),2) ).^2 ...
        + ( u(edges(:,1),2) - u(edges(:,2),2) ).^2;   
end

function result = integrate_u_ti_tj(gradUTotal, allSurroundingTriangles,...
    m, n)
%integrate_u_ti_tj berechnet int |nabla u|^2 T_i T_j.
    
    % komplette Matrix integrate_u_ti_tj
    result = zeros((m+1)*(n+1)); 
    
    %Diagonale der Matrix also die Einträge T_i T_i
    diag =zeros((m+1)*(n+1),1); 
    
    % erste Nebendiagonale von u, also die Einträge zu T_i T_i+1
    secDiag = zeros((m+1)*(n+1),1); 
   
    % n+1 te Nebendiagonale von u, die Einträge zu T_i T_i+n+1
    farDiag1 = zeros((m+1)*(n+1),1);     
    
    %n+2 te Nebendiagonale von u, die Einträge zu T_i, T_i+n+2
    farDiag2 = zeros((m+1)*(n+1),1); 
    
    for i=1:(m+1)*(n+1)
        % alle umliegenden Dreiecke zum Gitterpunkt i
        surroundings = allSurroundingTriangles(i,:); 
        % gibt gradU von den umliegenden Dreiecken an. 
        gradUSurroundings = zeros(6,1); 
        for j=1:6
            if (surroundings(j)~=0)
                gradUSurroundings(j)=gradUTotal(surroundings(j));
            end
        end
        
        % berechnet die diagonalen. Formeln wieder in der Bachelorarbeit. 
        diag(i)=(1/12) * sum(gradUSurroundings);  
        secDiag(i) = (1/24) * (gradUSurroundings(3) ...
            + gradUSurroundings(6) ) ; 
        farDiag1(i) = (1/24) * (gradUSurroundings(4) ...
            + gradUSurroundings(5) ) ;
        farDiag2(i) = (1/24) * (gradUSurroundings(5) ...
            + gradUSurroundings(6) ) ;
    end
    
    % wie vorhin ist die gewünschte Matrix symmetrisch und die Diagonalen
    % müssen auf beiden Seiten stehen. 
    secDiag2 =  [0; secDiag(1:(n+1)*(m+1)-1)]; 
    farDiag12 = [zeros(n+1,1) ; farDiag1(1:(n+1)*(m+1)-(n+1))] ;
    farDiag22 = [zeros(n+2,1) ; farDiag2(1:(n+1)*(m+1)-(n+2))] ;
    % schreibt alle Diagonalen an der richtigen Stelle in eine Sparse
    % Matrix
    result = spdiags([farDiag2 farDiag1 secDiag diag secDiag2 farDiag12 ...
        farDiag22 ],[-n-2 -n-1 -1 0 1 n+1 n+2],(n+1)*(m+1),(n+1)*(m+1) );    
end

function result = integrate_gti_gtj( allSurroundingTriangles, m, n)
%integrate_gti_gtj berechnet \int_\Omega nabla Ti nabla Tj dx
% Diese Matrix heißt später B
    
    % komplette Matrix integrate_u_ti_tj
    result = zeros((m+1)*(n+1)); 
    
    %Diagonale der Matrix
    diag =zeros((m+1)*(n+1),1); 
    
    % erste Nebendiagonale von u
    secDiag = zeros((m+1)*(n+1),1); 
   
    % n+1 te Nebendiagonale von u
    farDiag = zeros((m+1)*(n+1),1);     

    for i=1:(m+1)*(n+1)
        % berechnet alle umliegenden Dreiecke des Gitterpunktes i
        sur = allSurroundingTriangles(i,:);
        % alle Einträge werden 1 gesetzt, die vorher einen Wert hatten, da
        % für die Berechnung der Diagonalen für jeden Wert nur 1 eingesetzt
        % wird. So werden die nicht vorhandenen Dreiecke ausgelassen und
        % die anderen Berechnet. 
        for j=1:6
            if (sur(j)~=0) 
                sur(j)=1; 
            end
        end
        % Diagonalen werden berechnet. Formel wieder in der Bachelorarbeit
        diag(i)=(1/2) * sur(1) + (1/2) * sur(2) + sur(3) + sur(4) ...
            + (1/2) * sur(5) + (1/2) * sur(6);  
        secDiag(i) = -(1/2) * sur(3) - (1/2) * sur(6); 
        farDiag(i) =  -(1/2) * sur(4) - (1/2) * sur(5); 
    end
    % Da Matrix symmetrisch ist, müssen die Diagonalen auf beiden Seiten
    % der Matrix stehen. 
    secDiag2 =  [0; secDiag(1:(n+1)*(m+1)-1)]; 
    farDiag2 = [zeros(n+1,1) ; farDiag(1:(n+1)*(m+1)-(n+1))] ;
    % Ergebnisse werden in der Sparsematrix geschrieben. 
    result = spdiags([farDiag secDiag diag secDiag2 farDiag2 ]...
        ,[-n-1 -1 0 1 n+1],(n+1)*(m+1),(n+1)*(m+1) ); 
end

function result = integrate_ti_tj(allSurroundingTriangles, m, n)
%integrate_ti_tj berechnet integral T_i T_j dx
% Diese Matrix heißt D in der Ausarbeitung
    
    % komplette Matrix integrate_u_ti_tj
    result = zeros((m+1)*(n+1)); 
    
    %Diagonale der Matrix
    diag =zeros((m+1)*(n+1),1); 
    
    % erste Nebendiagonale von u
    secDiag = zeros((m+1)*(n+1),1); 
   
    % n+1 te Nebendiagonale von u
    farDiag1 = zeros((m+1)*(n+1),1);     
    
    % n+2 te Nebendiagonale von u
    farDiag2 = zeros((m+1)*(n+1),1); 
    
    for i=1:(m+1)*(n+1)
        % berechnet die umliegenden Dreiecke des Gitterpunktes i
        surroundings = allSurroundingTriangles(i,:); 
        % alle Einträge werden 1 gesetzt, die vorher einen Wert hatten, da
        % für die Berechnung der Diagonalen für jeden Wert nur 1 eingesetzt
        % wird. So werden die nicht vorhandenen Dreiecke ausgelassen und
        % die anderen Berechnet.       
        for j=1:6
            if (surroundings(j)~=0) 
                surroundings(j)=1; 
            end
        end
        % berechnet die Diagonalen. Formel in der Bachelorarbeit. 
        diag(i)=(1/12) * sum(surroundings);  
        secDiag(i) = (1/24) * surroundings(3) + (1/24) * surroundings(6); 
        farDiag1(i) =  (1/24) * surroundings(4) + (1/24) * surroundings(5); 
        farDiag2(i) =  (1/24) * surroundings(5) + (1/24) * surroundings(6); 
    end
    % Fast alle Diagonalen müssen zweimal eingetragen werden. 
    secDiag2 =  [0; secDiag(1:(n+1)*(m+1)-1)]; 
    farDiag12 = [zeros(n+1,1) ; farDiag1(1:(n+1)*(m+1)-(n+1))] ;
    farDiag22 = [zeros(n+2,1) ; farDiag2(1:(n+1)*(m+1)-(n+2))] ;
    % Diagonalen werden in der Sparsematrix geschrieben. 
    result = spdiags([farDiag2 farDiag1 secDiag diag secDiag2 farDiag12 ...
        farDiag22 ],[-n-2 -n-1 -1 0 1 n+1 n+2],(n+1)*(m+1),(n+1)*(m+1) );  
end

