function [edges, E] = triangulasation(m,n)
% TRIANGULASATION gibt die Triangulierung von dem Gebiet Omega an. 
%   edges ist eine n*m*2 x 3 Matrix und berechnet alle Ecken der Dreiecke
%   allSurroundingTriangles ist eine (m+1)*(n+1) x 6 Matrix und beinhaltet 
%   in der i ten Zeile alle indizes der Dreiecke, die um den Punkt i 
%   liegen. 
    
    edges = edges_of_triangles(m,n); 
    E = surrounding_triangles(m,n); 
end

function edges = edges_of_triangles(m,n)
% EDGES_OF_TRIANGLES berechnet die Ecken jedes Dreiecks. 

%   Nun möchte ich von meinem Gitter die Eckpunkt des jeweiligen Dreiecks 
%   bestimmen. 
%   _ _ _ _ _ _ 
%   |\ 2|\ 4|\ 6|
%   |1\ |3\ |5\ | 
%   |__\|__\|__\|
%   |\ 8|\10|\12|
%   |7\ |9\ |1\ | 
%   |__\|__\|__\|
%
%   hier gilt n=3, m=2
%   Zunächst berechnet man für alle Dreiecke den linken obersten Randpunkt.
%   Nun unterscheidet man zwischen geraden und ungeraden Dreiecken. 
%   Gerade Dreiecke: 
%   Die nächste Ecke liegt einfach rechts neben der schon berechneten Ecke, 
%   d.h. man rechnet den Wert +1
%   Die untere Ecke liegt direkt unter der vorherigen Ecke, also +n+1
%   Ungerade Dreiecke: 
%   Die nächste Ecke liegt einfach unter der schon berechneten Ecke, d.h. 
%   man rechnet den Wert +n+1
%   Die untere Ecke liegt direkt neben der vorherigen Ecke, also +1

edges=zeros(n*m*2,3);

for i=0:n*m*2-1
    edges(i+1,1)=floor(i/2)+floor(i/(2*n))+1;
    if (mod(i,2)~=0)
        edges(i+1,2)=edges(i+1,1)+1;
        edges(i+1,3)=edges(i+1,2)+n+1;
    else
        edges(i+1,2)=edges(i+1,1)+n+1;
        edges(i+1,3)=edges(i+1,2)+1;
    end
end

end


function E = surrounding_triangles(m,n)
%SURROUNDING_TRIANGLES gibt alle anliegenden Dreiecke von jedem Gitterpunkt 
%an.

%   Insgesamt brauchen wir die 6 Dreiecke, die um den Gitterpunkt j liegen.
%   Sobald wir das erste Dreieck (oben links) gefunden haben, können wir 
%   die anderen daraus berechnen. 
%
%   Berechnung der Lage des ersten Dreieckes:
%   Diese Formel scheint aus dem Himmel zu fallen, funktioniert aber. 

% initiierung von E
E = zeros((n+1)*(m+1),6);

for j=1:(n+1)*(m+1)
    if (j>n+1 && j<=(n+1)*m && mod(j,n+1)~=0 && mod(j,n+1)~=1)
        % Gitterpunkte in der Mitte
        E(j,1) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1));
        E(j,2)=E(j,1)+1;
        E(j,3)=E(j,1)+2;
        E(j,4)=E(j,1)+2*n+1;
        E(j,5)=E(j,4)+1;
        E(j,6)=E(j,4)+2;
    elseif (j==1) 
        % Gitterpunkt in der oberen linken Ecke
        E(j,5)=1; 
        E(j,6)=2;  
    elseif (j<n+1 && j>1)
        % oberer Rand
        E(j,4)= 2*j-(2*n+5)-2*floor((j-n-2)/(n+1)) + 2*n + 1;
        E(j,5)=E(j,4)+1;
        E(j,6)=E(j,4)+2;   
    elseif (j==n+1)
        % Gitterpunkt in der unteren linken Ecke
        E(j,4)= 2*n; 
    elseif (mod(j,n+1)==0 && j~= (n+1)*(m+1))
        % Gitterpunkte am rechten Rand
        E(j,1) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1));
        E(j,2)=E(j,1)+1;
        E(j,4)=E(j,1)+2*n+1;  
    elseif (j== (n+1)*(m+1))
        % Gitterpunkt an der rechten unteren Ecke
        E(j,1) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1));
        E(j,2)=E(j,1)+1;      
    elseif( mod(j,n+1)==1 && j~= (n+1)*m+1)
        % Gitterpunkte am linken Rand
        E(j,3)= 2*j-(2*n+5)-2*floor((j-n-2)/(n+1))+2;
        E(j,5)=E(j,3)+2*n;
        E(j,6)=E(j,5)+1;  
    elseif (j==(n+1)*m+1)
        % Gitterpunkt in der rechten oberen Ecke
        E(j,3) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1)) +2;           
    elseif (j>(n+1)*m+1 && j< (n+1)*(m+1)) 
        % unterer Rand
        E(j,1) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1));
        E(j,2)=E(j,1)+1;
        E(j,3)=E(j,1)+2;        
    end
end


end
