function result = integrate_u_ti_tj(gradUTotal, allSurroundingTriangles, m, n)
%integrate_u_ti_tj berechnet int |nabla u|^2 T_i T_j. Diese Matrix heiﬂt A
%in der Ausarbeitung
    
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
        gradUSurroundings = grad_u_surroundings(surroundings, gradUTotal);
        diag(i)=(1/12) * sum(gradUSurroundings);  
        secDiag(i) = (1/24) * (gradUSurroundings(3) + gradUSurroundings(6) ) ; 
        farDiag1(i) = (1/24) * (gradUSurroundings(4) + gradUSurroundings(5) ) ;
        farDiag2(i) = (1/24) * (gradUSurroundings(5) + gradUSurroundings(6) ) ;
    end
    
    secDiag2 =  [0; secDiag(1:(n+1)*(m+1)-1)]; 
    farDiag12 = [zeros(n+1,1) ; farDiag1(1:(n+1)*(m+1)-(n+1))] ;
    farDiag22 = [zeros(n+2,1) ; farDiag2(1:(n+1)*(m+1)-(n+2))] ;
    
    result = spdiags([farDiag2 farDiag1 secDiag diag secDiag2 farDiag12 farDiag22 ],[-n-2 -n-1 -1 0 1 n+1 n+2],(n+1)*(m+1),(n+1)*(m+1) );  
    
    
    
end

