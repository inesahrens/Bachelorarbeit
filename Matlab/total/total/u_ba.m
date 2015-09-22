n=6;
m=6;
nu = .010; 
epsilon = [0.01;nu* 20*1/n ;(20*1/n )/nu]; 

v = ones((m+1)*(n+1),1); 
for i=0:m
    v(i*(n+1)+n/2) = 0;
    v(i*(n+1)+n/2+1) = 0;
end

u0=zeros((n+1)*(m+1), 1);
for i=0:m
     u0(i*(n+1)+1) = 1;
     u0(i*(n+1)+n+1) = 2;
end
u0(1)= 2;
u0(m+1)=3;
u0(m+2)=2;
u0(2*m+2)=3;
u0((m+1)*(n-1)+1)=2
u0((m+1)*n)=3
u0((m+1)*n+1)=2
u0((m+1)*(n+1))=3
a = reshape(u0,n+1,m+1)

[edges, allSurroundingTriangles] = triangulasation(m,n); 



u = calculate_u(edges, allSurroundingTriangles, u0, epsilon(1), v, m, n); 

u = reshape(u,n+1,m+1); 

[X,Y] = meshgrid(0:1/n:1);
figure 
mesh(X,Y,u)