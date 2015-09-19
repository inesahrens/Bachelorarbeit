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

