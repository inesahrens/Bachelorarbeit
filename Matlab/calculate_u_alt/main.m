%Deklaration der Variablen: 

% n+1 ist die Anzahl der Stützstellen nach rechts
% m+1 ist die Anzahl der Stützstellen nach unten
n=4;
m=4;

epsilon =0.01;
% v ist die Auswertung der Funktion v an den Stützstellen 1,...,(n+1)(m+1)
% später wird v berechnet. Noch ist es einfach eine festgelegte Funktion
v = ones((m+1)*(n+1), 1); 
for i=0:m
    v(i*(n+1)+n/2+1) = 0;
end

% u0 ist die erweiterte Funktion der Randbedingung. Sie ist auf dem Rand
% Gamma 1 und Gamma 2 u0 und sonst 0. 
% sie wird hier auch erstmal nur gesetzt, später gibt es andere Werte. 
u0=sparse((n+1)*(m+1), 1);
for i=0:m
    u0(i*(n+1)+1) = 2;
    u0(i*(n+1)+n+1) = 2;
end

% Für die Triangulierung bestimme ich zunächst alle Eckpunkte der Dreiecke
% und alle benachbarten Dreiecke von einem Gitterpunkt.
% edgesOfTriangles ist eine n*m*2 x 3 Matrix 
% surroundingTriangles ist ein (m+1)(n+1) x 6 Matrix  
[edges, allSurroundingTriangles] = triangulasation(m,n); 

% Da für die Berechnung des Integrals für alle Dreiecke das Integral von
% v^2+eps gebraucht wird, berechnen wir es hier. 
integral_v = integral_of_v_total(v, edges, epsilon,m,n);

% Ziel ist es  
%  - \int_Omega (v^2 + eps) \nabla u0 \nabla T_j 
% = \sum_i^{(n+1)(m+1)} \int_Omega (v^2+eps) \nabla T_i \nabla T_j
% für alle j =1,... (m+1)(n+1) zu berechnen

% linke Seite:
% die linke Seite ist 0 für alle Werte im Inneren, bis auf die Werte, die
% direkt am Rand sind. 
% die Werte am Rand und am 2. Rand müssen berechnet werden. 
 
links = - integrate_u0(edges, allSurroundingTriangles, integral_v, u0, m,n ); 

% rechte Seite
% die rechte Seite ist 0 für alle Werte am Rand und sonst die Werte des
% Integrals. 

rechts = integrate_u(allSurroundingTriangles,integral_v, n, m); 
% hier wird x aus  Ax=b berechnet. 
result = rechts\links; 

result = reshape(result,n+1,m+1); 
u0 = reshape(u0,n+1,m+1); 
result = result + u0; 


% Plotten des Ergebnisses 
% noch geht das nur für Gitter, die gleichhoch wie breit sind. 
[X,Y] = meshgrid(0:1/n:1);
figure
mesh(X,Y,result)

