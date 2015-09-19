function edges = edges_of_triangles(m,n)
%integrate_of_v_square berechnet \int_E_i v^2 
% dabei gilt:
% v ist der vektor, der den Riss darstellt 
% m ist die Zeilenanzahl des triangulierten Gebietes
% n ist die Spaltenanzahl des triangulierten Gebietes
% E_i das i-te Dreieck

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

