function integrate = integral_of_v_total(v, edges, epsilon,m,n)
%integral_of_v_total berechnet int v^2 + epsilon über alle dreiecke
integrate = zeros(n*m*2,1);
for i=1:n*m*2
    integrate(i) = integral_of_v(v(edges(i,:)),epsilon);
end

end

function integral  = integral_of_v(v, epsilon)
%integrate_of_v berechnet \int_E_i v^2 + epsilon
% dabei gilt:
% eps ist der kleine Parameter Epsilon
% vi ist der Wert von v an einem Eckpunkt
% Die Formel für das Integral wurde in der zugehörigen Bachelorarbeit
integral= 1/12*(v(1)^2 + v(2)^2 + v(3)^2 + v(1)*v(2) + v(1)*v(3) + v(2)*v(3)) + 1/2 * epsilon; 

end
