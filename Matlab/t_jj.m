function Tjj = t_jj(v,m,n,j,epsilon)
% Berechnung von T_jj (bessere Bezeichnung!) 

% Berechnung der Indizes düe die Eckpunkte aller Dreiecke

[index_v1, index_v2,index_v3]= edges_of_triangles(m,n);

% Berechnung von dem Integral über das Dreieck E_i

% Insgesamt brauchen wir die 6 Dreiecke, die um den Gitterpunkt j liegen.
% Sobald wir das erste Dreieck (wieder oben links) gefunden haben, können wir die anderen
% daraus berechnen. 
%
% Berechnung die umliegenden 6 Dreiecke
[E1,E2,E3,E4,E5,E6]=Surrounding_Triangles(n,j); 

% Berechnung des integrals über das Dreieck i. Dabei ist v(index_v1(i)) der Wert von v an der linken oberen Ecke des Dreiecks i.
integral_dreieck_1 = integrate_of_v_square(v(index_v1(E1)), v(index_v2(E1)), v(index_v3(E1)), epsilon);
integral_dreieck_2 = integrate_of_v_square(v(index_v1(E2)), v(index_v2(E2)), v(index_v3(E2)), epsilon);
integral_dreieck_3 = integrate_of_v_square(v(index_v1(E3)), v(index_v2(E3)), v(index_v3(E3)), epsilon);
integral_dreieck_4 = integrate_of_v_square(v(index_v1(E4)), v(index_v2(E4)), v(index_v3(E4)), epsilon);
integral_dreieck_5 = integrate_of_v_square(v(index_v1(E5)), v(index_v2(E5)), v(index_v3(E5)), epsilon);
integral_dreieck_6 = integrate_of_v_square(v(index_v1(E6)), v(index_v2(E6)), v(index_v3(E6)), epsilon);

% Berechnung von T_jj
% Die Formel dazu ist in der zugehörigen Bachelorarbeit
% Darf nur ausgerechet werden, falls wir nicht am Rand sind. Für den Rand müssen andere Dinge gelten! 
Tjj = 2*integral_dreieck_3+ 1*integral_dreieck_6 + integral_dreieck_5 + 2*integral_dreieck_4+1*integral_dreieck_1 + 1*integral_dreieck_2;


