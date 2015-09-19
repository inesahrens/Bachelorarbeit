function result = integrate_ti_tj(allSurroundingTriangles, m, n)
%integrate_ti_tj berechnet integral T_i T_j dx
% Diese Matrix heiﬂt A in der Ausarbeitung
    
    % komplette Matrix integrate_u_ti_tj
    result = zeros((m+1)*(n+1)); 
    
    %Diagonale der Matrix
    diag =zeros((m+1)*(n+1),1); 
    
    % berechnet die erste Nebendiagonale von u
    secDiag = zeros((m+1)*(n+1),1); 
   
    % berechnet die n+1 te Nebendiagonale von u
    farDiag1 = zeros((m+1)*(n+1),1);     
    
    % berechnet die n+2 te Nebendiagonale von u
    farDiag2 = zeros((m+1)*(n+1),1); 
    
    for i=1:(m+1)*(n+1)
        surroundings = allSurroundingTriangles(i,:); 
        for j=1:6
            if (surroundings(j)~=0) 
                surroundings(j)=1; 
            end
        end
        diag(i)=(1/12) * sum(surroundings);  
        secDiag(i) = (1/24) * surroundings(3) + (1/24) * surroundings(6); 
        farDiag1(i) =  (1/24) * surroundings(4) + (1/24) * surroundings(5); 
        farDiag2(i) =  (1/24) * surroundings(5) + (1/24) * surroundings(6); 
    end
    
    secDiag2 =  [0; secDiag(1:(n+1)*(m+1)-1)]; 
    farDiag12 = [zeros(n+1,1) ; farDiag1(1:(n+1)*(m+1)-(n+1))] ;
    farDiag22 = [zeros(n+2,1) ; farDiag2(1:(n+1)*(m+1)-(n+2))] ;
    
    result = spdiags([farDiag2 farDiag1 secDiag diag secDiag2 farDiag12 farDiag22 ],[-n-2 -n-1 -1 0 1 n+1 n+2],(n+1)*(m+1),(n+1)*(m+1) );  
end

