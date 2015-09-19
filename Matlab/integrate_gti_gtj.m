function result = integrate_gti_gtj( allSurroundingTriangles, m, n)
%integrate_gti_gtj berechnet integral nabla Ti nabla Tj dx
% Diese Matrix heiﬂt B in der Ausarbeitung
    
    % komplette Matrix integrate_u_ti_tj
    result = zeros((m+1)*(n+1)); 
    
    %Diagonale der Matrix
    diag =zeros((m+1)*(n+1),1); 
    
    % berechnet die erste Nebendiagonale von u
    secDiag = zeros((m+1)*(n+1),1); 
   
    % berechnet die n+1 te Nebendiagonale von u
    farDiag = zeros((m+1)*(n+1),1);     

    
    for i=1:(m+1)*(n+1)
        sur = allSurroundingTriangles(i,:);
        for j=1:6
            if (sur(j)~=0) 
                sur(j)=1; 
            end
        end
        diag(i)=(1/2) * sur(1) + (1/2) * sur(2) + sur(3) + sur(4) + (1/2) * sur(5) + (1/2) * sur(6);  
        secDiag(i) = -(1/2) * sur(3) - (1/2) * sur(6); 
        farDiag(i) =  -(1/2) * sur(4) - (1/2) * sur(5); 
    end
    secDiag2 =  [0; secDiag(1:(n+1)*(m+1)-1)]; 
    farDiag2 = [zeros(n+1,1) ; farDiag(1:(n+1)*(m+1)-(n+1))] ;
    
    result = spdiags([farDiag secDiag diag secDiag2 farDiag2 ],[-n-1 -1 0 1 n+1],(n+1)*(m+1),(n+1)*(m+1) );  

end

