function result = integrate_u(allSurroundingTriangles, integral_v, m, n)
%integrate_u berechnet das Integral über u, also  \int_Omega (v^2+eps)
%\nabla T_i \nabla T_j dx

    
    % result ist eine Sparse Matrix, die nur auf der Diagonalen, der ersten
    % Nebendiagonalen und der n+1 ten Nebendiagonalen Einträge hat.
    
    % berechnet die Diagonale von u. Hier werden Einträge der Form 
    % \int_Omega (v^2+eps) \nabla T_i \nabla T_i dx berechnet
    diag = diagonal_u(allSurroundingTriangles,integral_v, n, m); 
    % berechnet die erste Nebendiagonale von u. Also Einträge der Form 
    % \int_Omega (v^2+eps) \nabla T_i \nabla T_i+1 dx berechnet
    secDiag = second_diagonal_u(allSurroundingTriangles,integral_v, n, m);  
    % berechnet die n+1 te Nebendiagonale von u, also Einträge der Form
    % \int_Omega (v^2+eps) \nabla T_i \nabla T_i+1+n dx berechnet
    farDiag = far_diagonal_u(allSurroundingTriangles,integral_v, n, m); 
    % damit spDiags auf die richtigen Elemente zugreift, muss der erste bzw
    % die ersten n+1 Elemente 0 gesetzt werden. 
    secDiag2 =  [0; secDiag(1:(n+1)*(m+1)-1)];
    farDiag2 = [zeros(n+1,1) ; farDiag(1:(n+1)*(m+1)-(n+1))] ;
    % schreibt alle Diagonalen in eine Sparse Matrix
    result = spdiags([farDiag secDiag diag secDiag2 farDiag2 ]...
        ,[-n-1 -1 0 1 n+1],(n+1)*(m+1),(n+1)*(m+1) );  
end


function result = diagonal_u(allSurroundingTriangles, integral_v, n, m)
%DIAGONAL_U berechnet \int_Omega (v^2+eps) \nabla T_i \nabla T_i dx 

    result = zeros((n+1)*(m+1),1) ; 
    
    for i=1:(m+1)*(n+1)
        % berechnet die umliegenden Dreiecke vom Gitterpunkt i
        surrounders = allSurroundingTriangles(i,:); 
        % berechnet das \int_Omega (v^2+eps) dx für die umliegenden
        % Dreiecke
        intvsurround = integral_of_v_surroundings(surrounders, integral_v); 
        % rechnet die Integrale mit den richtigen Vorfaktoren zusammen. Die
        % Berechnung der Vorfaktoren kann in der Bachelorarbeit nachgelesen
        % werden. 
        result(i) =  intvsurround(1) + intvsurround(2) ...
            + 2 * intvsurround(3) + 2 * intvsurround(4) ...
            + intvsurround(5) + intvsurround(6) ;    
    end   
end


function result = second_diagonal_u(allSurroundingTriangles,integral_v, n, m)
%SECOND_DIAGONAL_U berechnet \int_Omega (v^2+eps) \nabla T_i \nabla T_i+1 dx 

    result = zeros((n+1)*(m+1),1); 
    
    for i=1:(n+1)*(m+1)
       % berechnet die umliegenden Dreiecke vom Gitterpunkt i
       surrounders = allSurroundingTriangles(i,:);
       % berechnet das \int_Omega (v^2+eps) dx für die umliegenden
       % Dreiecke
       intvsurround = integral_of_v_surroundings(surrounders, integral_v); 
       % rechnet die Integrale mit den richtigen Vorfaktoren zusammen. Die
       % Berechnung der Vorfaktoren kann in der Bachelorarbeit nachgelesen
       % werden.        
       result(i) = - intvsurround(3) - intvsurround(6);    
    end
end


function result = far_diagonal_u(allSurroundingTriangles,integral_v, n, m)
%FAR_DIAGONAL_U berechnet \int_Omega (v^2+eps) \nabla T_i \nabla T_i+1 dx 

    result = zeros((n+1)*(m+1),1); 
    
    for i=1:(n+1)*(m+1)
       % berechnet die umliegenden Dreiecke vom Gitterpunkt i
       surrounders = allSurroundingTriangles(i,:);
       % berechnet das \int_Omega (v^2+eps) dx für die umliegenden
       % Dreiecke       
       intvsurround = integral_of_v_surroundings(surrounders, integral_v);
       % rechnet die Integrale mit den richtigen Vorfaktoren zusammen. Die
       % Berechnung der Vorfaktoren kann in der Bachelorarbeit nachgelesen
       % werden.        
       result(i) = - intvsurround(4) - intvsurround(5);    
    end
end




