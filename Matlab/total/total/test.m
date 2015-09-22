%elements=[1,4,5;1,2,5;2,5,6;2,3,6;4,7,8;4,5,8;5,8,9;5,6,9];
%nodes = [1,1;1,2;1,3;2,1;2,2;2,3;3,1;3,2;3,3];

m=4;
n=4;
u=zeros((n+1)*(m+1),2);
[edges, allSurroundingTriangles] = triangulasation(m,n); 
[A, B, c, D] =  Matrixes_for_G(u, edges, allSurroundingTriangles, m, n ); 

nodes = zeros((n+1)*(m+1),2);
for i=1:(n+1)*(m+1)
    nodes(i,1)= floor(i/(n+1.001))+1; 
    nodes(i,2)= mod(i,(n+1)); 
    if nodes(i,2)==0
        nodes(i,2)=n+1;
    end
end
nodes

[M,L] = massStiffMatrix(nodes,edges);
full(L);
full(B);
if (B==L)
    b=1;
else
    b=0;
end

d=zeros((n+1)*(m+1));
for i=1:(n+1)*(m+1)
    for j=1:(n+1)*(m+1)
        if (M(i,j)==D(i,j))
            d(i,j)=1;
        end
    end
end
e = max(abs(M-D),0)
b
d
M(7,7)
D(7,7)
comp = full(M);
mein = full(D);
