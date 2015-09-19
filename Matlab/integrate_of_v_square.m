function [ integral ] = integrate_of_v_square(v, epsilon)
%integrate_of_v_square berechnet \int_E_i v^2 + epsilon
% dabei gilt:
% eps ist der kleine Parameter Epsilon
% vi ist der Wert von v an einem Eckpunkt
% Die Formel für das Integral wurde in der zugehörigen Bachelorarbeit
% berechnet

integral= 1/12*(v(1)^2 + v(2)^2 + v(3)^2 + v(1)*v(2) + v(1)*v(3) + v(2)*v(3)) + 1/2*epsilon;

end

