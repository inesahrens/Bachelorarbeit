function result = integrate_u(allSurroundingTriangles, integral_v, n, m)
%integrate_u berechnet das Integral über u, also \sum_i^{(n+1)(m+1)} \int_Omega (v^2+eps) \nabla T_i \nabla T_j

    
    % result ist eine Sparse Matrix, die nur auf der Diagonalen, der ersten
    % Nebendiagonalen und der n+1 ten Nebendiagonalen Einträge hat.
    
    % berechnet die Diagonale von u
    diag = diagonal_u(allSurroundingTriangles,integral_v, n, m); 
    % berechnet die erste Nebendiagonale von u
    secDiag = second_diagonal_u(allSurroundingTriangles,integral_v, n, m); 
    % damit spDiags auf die richtigen Elemente zugreift, muss der erste
    % wert 0 gesetzt werden. 
    secDiag2 =  [0; secDiag(1:(n+1)*(m+1)-1)]; 
    % berechnet die n+1 te Nebendiagonale von u
    farDiag = far_diagonal_u(allSurroundingTriangles,integral_v, n, m); 
    farDiag2 = [zeros(n+1,1) ; farDiag(1:(n+1)*(m+1)-(n+1))] ;
    result = spdiags([farDiag secDiag diag secDiag2 farDiag2 ],[-n-1 -1 0 1 n+1],(n+1)*(m+1),(n+1)*(m+1) );  
     
     % Da u=0 auf dem Rand ist, müssen alle Zeilen, die den Rand darstellen
     % 0 gesetzt werden, bis auf den Eintrag, in dem 
     
%      result(1:n+1:((n+1)*(m+1)),:)=0; 
%      result((1:n+1:((n+1)*(m+1)))+n,:)=0;
%      result(1:n+1:((n+1)*(m+1)), 1:n+1:((n+1)*(m+1)))=1; 
%      result((1:n+1:((n+1)*(m+1)))+n,(1:n+1:((n+1)*(m+1)))+n)=1;
      for i=0:m
         result(i*(n+1)+1,: ) = 0;
         result(i*(n+1)+1+n,:) = 0;        
         result(i*(n+1)+1,i*(n+1)+1) = 1;
         result(i*(n+1)+1+n,i*(n+1)+1+n ) = 1;
      end

end

