function gradU = grad_u_total(u, edges, m, n )
%grad_u_total berechnet |\nabla u|^2 für alle Dreiecke der Triangulierung. 
% Dabei ist  |\nabla u|^2 = (u31-u21)^2 + (u31-u21)^2 + (u32-u22)^2 +
% (u32-u22)^2
    gradU = ( u(edges(:,3),1) - u(edges(:,2),1) ).^2 + ( u(edges(:,1),1) - u(edges(:,2),1) ).^2 +  ( u(edges(:,3),2) - u(edges(:,2),2) ).^2 +  ( u(edges(:,1),2) - u(edges(:,2),2) ).^2;   

end

