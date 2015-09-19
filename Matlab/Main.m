
clc;
clear all; 

n=20;
m=20;
epsilon=1; 
%v=rand([(m+1)*(n+1) 1]);
v=ones([(m+1)*(n+1) 1]); 
u0=ones([(m+1)*(n+1) 1])-0.5; 

edges = edges_of_triangles(m,n);
% Bestimmt die Eckpunkte jedes Dreieckes. Dabei sind die Dreiecke wie folgt
% nummieriert: 
% _ _ _ _ _ _ 
%|\ 2|\ 4|\ 6|
%|1\ |3\ |5\ | 
%|__\|__\|__\|
%|\ 8|\10|\12|
%|7\ |9\ |1\ | 
%|__\|__\|__\|
%

integrate = zeros(n*m*2+1,1);
for i=1:n*m*2
    integrate(i+1) = integrate_of_v_square(v(edges(i,:)),epsilon);
end
 
% bestimmt alle Integrale über alle Dreiecke

E = Surrounding_Triangles(m,n);
% Bestimmt alle 6 Dreiecke, die um den jeden inneren Punkt liegen. Für
% Randpkt ist der Term 0. 

A= zeros((n+1)*(m+1),(n+1)*(m+1)) ;
B = zeros((n+1)*(m+1),(n+1)*(m+1)) ;
% die Schleife sortiert die Randwerte aus. Im If Teil stehen alle Werte im
% Inneren und im Else Teil stehen alle Randwerte. 
for j=1:(n+1)*(m+1)
        % Hauptdiagonale
        A(j,j)= 2*integrate(E(j,3)+1) + integrate(E(j,6)+1) +  integrate(E(j,5)+1) + 2*integrate(E(j,4)+1) + integrate(E(j,1)+1) + integrate(E(j,2)+1);
        B(j,j)=A(j,j); 
        % erste Nebendiagonale
        if (j<(n+1)*(m+1))
            A(j,j+1)=-1*integrate(E(j,3)+1) + -1*integrate(E(j,6)+1);
            A(j+1,j)= A(j,j+1); 
        end
        % Nebendiagonale weit weg
        if (j<(n)*(m+1))
            A(j,j+1+n)= -1*integrate(E(j,4)+1) + -1*integrate(E(j,5)+1);
            A(j+1+n,j)= A(j,j+1+n); 
        end
end

u= rand([(m+1)*(n+1) 1]);


for k=1:10000
   u=u0_vektor(m,n,u,u0);
   s=transpose ( (A*u)\B) ; 
   u=s + u; 
   c(:,k)=u; 
end







