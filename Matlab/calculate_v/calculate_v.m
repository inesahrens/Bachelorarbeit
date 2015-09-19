function result =calculate_v( )
%calculate_v berechnet v. 
% das verfahren wird ein Newtonverfahren

% n+1 ist die Anzahl der Stützstellen nach rechts
% m+1 ist die Anzahl der Stützstellen nach unten

n=20;
m=20;

epsilon = 0.1;

v = .5*ones((m+1)*(n+1),1); 
for i=0:m
    v(i*(n+1)+n/2) = 0;
end

eta = ones((n+1)*(m+1),1);  
v0 = ones((n+1)*(m+1),1); 
const = 1; 
%hier wird fürs erste u gesetzt. Später bekommt diese Funktion u gegeben
u=ones((n+1)*(m+1),2);
%u(:,1)=[1;2;3;2;3;4;3;4;5];
%u(:,2)=[1;2;3;2;3;4;3;4;5];
%u(:,1)=[1;2;3;4;5;2;3;4;5;6;3;4;5;6;7;4;5;6;7;8;5;6;7;8;9];
%u(:,2)=[1;2;3;4;5;2;3;4;5;6;3;4;5;6;7;4;5;6;7;8;5;6;7;8;9];

% edges ist eine n*m*2 x 3 Matrix und berechnet alle Ecken der Dreiecke
% allSurroundingTriangles ist eine (m+1)*(n+1) x 6 Matrix Beinhaltet in der
% i ten Zeile alle indizes der Dreiecke, die um den Punkt i liegen. 
[edges, allSurroundingTriangles] = triangulasation(m,n); 

%Berechnet alle Matrizen, die für die berechnung von G gebraucht werden,
%die unabhängig von v sind. 
[A, B, c, D] = Matrixes_for_G(u, edges, allSurroundingTriangles, m, n ); 

%berechnung von G1v
G1v = 2*A + 2*epsilon*B + (2/epsilon)*D; 

%berechnung von G1eta
G1eta = D; 

for i=1:50
%while (norm(v1-v)>.00002)

    %berechnet das Komplette Integral
    %todo intEta lokal in intetatitotal schreiben. 
    intEtaTotal = integrate_eta_ti_total(edges, allSurroundingTriangles, eta, m, n ); 
    
    % berechnung von G1
    G1 = (2*A+2*epsilon*B+ (2/epsilon)*D)*v - (2/epsilon)*c + intEtaTotal; 
    
    [G2, G2v, G2eta] = G2_and_deratives(v, eta, const, v0, D, edges, allSurroundingTriangles, m, n); 
    gradG = sparse ((m+1)*(n+1)*2 ,(m+1)*(n+1)*2 ) ; 

    gradG(1:(m+1)*(n+1),1:(m+1)*(n+1)) = G1v; 
    gradG(1:(m+1)*(n+1),(m+1)*(n+1)+1:(m+1)*(n+1)*2) = G1eta; 
    gradG((m+1)*(n+1)+1:(m+1)*(n+1)*2,1:(m+1)*(n+1)) = G2v;     
    gradG((m+1)*(n+1)+1:(m+1)*(n+1)*2,(m+1)*(n+1)+1:(m+1)*(n+1)*2) = G2eta; 
    G = sparse ((m+1)*(n+1)*2,1) ;
    G(1:(m+1)*(n+1),1) = - G1;
    G((m+1)*(n+1)+1:(m+1)*(n+1)*2,1) = - G2; 
    
    %[G,gradG] = calculate_G_gradG(edges, allSurroundingTriangles, v, eta, const, v0, m, n  )
    
    s = gradG\G ; 
    v = s(1:(m+1)*(n+1)) + v; 
    eta = s((m+1)*(n+1)+ 1:(m+1)*(n+1)*2) + eta; 
%     resultv = reshape(v,n+1,m+1);
%     [X,Y] = meshgrid(0:1/n:1);
% figure
% mesh(X,Y,resultv)
end

resultv = reshape(v,n+1,m+1); 

resulteta = reshape(eta,n+1,m+1);

% Plotten des Ergebnisses 
% noch geht das nur für Gitter, die gleichhoch wie breit sind. 
[X,Y] = meshgrid(0:1/n:1);
figure
mesh(X,Y,resultv)

figure 
mesh(X,Y,resulteta)
end

