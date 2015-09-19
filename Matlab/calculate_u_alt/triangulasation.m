function [edges, E] = triangulasation(m,n)
    
    edges = edges_of_triangles(m,n); 
    E = surrounding_triangles(m,n); 

end



function edges = edges_of_triangles(m,n)

%Nun möchte ich von meinem Gitter die Eckpunkt des jeweiligen Dreiecks bestimmen. 
% _ _ _ _ _ _ 
%|\ 2|\ 4|\ 6|
%|1\ |3\ |5\ | 
%|__\|__\|__\|
%|\ 8|\10|\12|
%|7\ |9\ |1\ | 
%|__\|__\|__\|
%
% hier gilt n=3, m=2
% Zunächst berechnet man für alle Dreiecke den linken obersten Randpunkt.
% Nun unterscheidet man zwischen geraden und ungeraden Dreiecken. 
% Gerade Dreiecke: 
% Die nächste Ecke liegt einfach rechts neben der schon berechneten Ecke, d.h. man rechnet den Wert +1
% Die untere Ecke liegt direkt unter der vorherigen Ecke, also +n+1
% Ungerade Dreiecke: 
% Die nächste Ecke liegt einfach unter der schon berechneten Ecke, d.h. man
% rechnet den Wert +n+1
% Die untere Ecke liegt direkt neben der vorherigen Ecke, also +1

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
%Surrounded_Triangles gibt alle Integrale über alle anliegenden Dreiecke

% Berechnung der Indizes düe die Eckpunkte aller Dreiecke

%%[index_v1, index_v2,index_v3]= edges_of_triangles(m,n);

% Insgesamt brauchen wir die 6 Dreiecke, die um den Gitterpunkt j liegen.
% Sobald wir das erste Dreieck (wieder oben links) gefunden haben, können wir die anderen
% daraus berechnen. 
%
% Berechnung der Lage des ersten Dreieckes:
% Diese Formel scheint aus dem Himmel zu fallen, funktioniert aber. 
E=zeros((n+1)*(m+1),6);

for j=1:(n+1)*(m+1)
    if (j>n+1 && j<=(n+1)*m && mod(j,n+1)~=0 && mod(j,n+1)~=1)
        E(j,1) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1));
        E(j,2)=E(j,1)+1;
        E(j,3)=E(j,1)+2;
        E(j,4)=E(j,1)+2*n+1;
        E(j,5)=E(j,4)+1;
        E(j,6)=E(j,4)+2;
    elseif (j==1) 
        E(j,5)=1; 
        E(j,6)=2;  
    elseif (j<n+1 && j>1)
        E(j,4)= 2*j-(2*n+5)-2*floor((j-n-2)/(n+1)) + 2*n + 1;
        E(j,5)=E(j,4)+1;
        E(j,6)=E(j,4)+2;   
    elseif (j==n+1)
        E(j,4)= 2*n; 
    elseif (mod(j,n+1)==0 && j~= (n+1)*(m+1))
        E(j,1) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1));
        E(j,2)=E(j,1)+1;
        E(j,4)=E(j,1)+2*n+1;  
    elseif (j== (n+1)*(m+1))
        E(j,1) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1));
        E(j,2)=E(j,1)+1;      
    elseif( mod(j,n+1)==1 && j~= (n+1)*m+1)
        E(j,3)= 2*j-(2*n+5)-2*floor((j-n-2)/(n+1))+2;
        E(j,5)=E(j,3)+2*n;
        E(j,6)=E(j,5)+1;  
    elseif (j==(n+1)*m+1)
        E(j,3) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1)) +2;           
    elseif (j>(n+1)*m+1 && j< (n+1)*(m+1)) 
        E(j,1) = 2*j-(2*n+5)-2*floor((j-n-2)/(n+1));
        E(j,2)=E(j,1)+1;
        E(j,3)=E(j,1)+2;        
    end
end


end
