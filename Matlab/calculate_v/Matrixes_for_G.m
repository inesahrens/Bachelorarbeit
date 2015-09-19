function [ G2, G2v, G2eta ] = G2_and_deratives(v, eta, const, v0, D, edges, allSurroundingTriangles, m, n)
%G2_AND_DERATIVES berechnet G2 und die Ableitung nach v und eta. 
%   G2(v, eta) = \int_\Omega f(v,eta) T_i dx \forall i 
%   mit f(v,eta) = eta - max{0,eta + c(v-v0)} - min{0,eta+cv} 
%   oder anders geschrieben
%   f(v,eta) = -c(v-v0) falls -c(v-v0) <= eta
%              eta      falls -cv < eta <  -c(v-v0) 
%              -cv      falls  eta <=  -cv 
% 
%    f_v(v,eta) =  -c       falls -c(v-v0) < eta oder eta < -cv 
%                  0        falls -cv < eta <  -c(v-v0) 
%                 [-c,0]    falls -c(v-v0) = eta oder eta = -cv 
% 
%    f_eta(v,eta)= 0        falls -c(v-v0) < eta oder eta < -cv 
%                  1        falls -cv < eta <  -c(v-v0) 
%                 [0,1]     falls -c(v-v0) = eta oder eta = -cv 
%   damit ist die gesamte Ableitung
%   G_2v(v, eta)   = \int_\Omega f_v   T_i T_j dx
%   G_2eta(v, eta) = \int_\Omega f_eta T_i T_j dx
    
    % Berechnet G2
    G2 = calculate_G2(v, eta, const, v0,D, m, n); 
    % Berechnet die ableitungen von f
    [f_v, f_eta] = calculate_f_v_eta(v, eta, const, v0, m, n); 
    % berechnet die gesamten Ableitungen
    G2v = calculate_G2div(f_v, edges, allSurroundingTriangles, m, n); 
    G2eta = calculate_G2div(f_eta, edges, allSurroundingTriangles, m, n); 
end
  
function G2 = calculate_G2(v, eta, const, v0, D, m, n)

    % initialisierung
    w = zeros((n+1)*(m+1),1); 
    % alle Fallunterscheidungen in der Formel sind realisiert
    for i=1:(n+1)*(m+1)
        if ( -const*(v(i)-v0(i)) <= eta(i) ) 
            w(i) = -const*(v(i)-v0(i)) ; 
        elseif ( - const * v(i) < eta(i) && eta(i) <  -const*(v(i)-v0(i)) )
            w(i) = eta(i); 
        elseif ( eta(i) <= -const* v(i) )
            w(i) = - const * v(i) ; 
        end
    end
    % alle Vorfaktoren aus dem min und max können aus dem integral gezogen
    % werden. Damit bleibt nur die Matrix D = \int_\Omega T_i T_j zu berechnen,
    % die vorher schon berechnet wurde. 
    G2 = D * w; 

end

function [f_v, f_eta] = calculate_f_v_eta(v, eta, const, v0, m, n)

    f_v = zeros((n+1)*(m+1),1); 
    f_eta = zeros((n+1)*(m+1),1); 
    % hier werden einfach wiede die Fallunterscheidungen aufgenommen, die
    % oben bereits erklärt wurden. 
    for i=1:(n+1)*(m+1)
        for j=i:(n+1)*(m+1)
            if ( -const*(v(i)-v0(i)) < eta(i) || eta(i) < -const* v(i)) 
                f_v(i) = -const ; 
                f_eta(i) = 0 ; 
            elseif ( - const * v(i) < eta(i) && eta(i) <  -const*(v(i)-v0(i)) )
                f_v(i) = 0; 
                f_eta(j) = 1; 
            elseif ( -const*(v(i)-v0(i)) == eta(i) || eta(i) == -const* v(i) )
                  % es f_v und f_eta aus dem intervall gewählt werden. Der
                  % Einfachheit halber habe ich jeweils 0 gewählt. 
                  f_v(i) = 0; 
                  f_eta(i) = 0; 
            end
        end
    end  
end

function G2div = calculate_G2div(f, edges, allSurroundingTriangles, m, n)
% CALCULATE_G2DIV berechnet die Ableitung von G_2. 
%   Es ist egal, ob die ableitung nach v oder eta betrachet wird, da die unterschiede nur in f liegen. f wurde bereits berechnet.   

    % initialisierung
    diag = zeros((n+1)*(m+1),1); 
    secdiag = zeros((n+1)*(m+1),1); 
    thirddiag = zeros((n+1)*(m+1),1);
    forthdiag = zeros((n+1)*(m+1),1); 
    for i=1:(n+1)*(m+1)
       % umliegende Dreieck des Gitterpunktes i berechnen 
       surrounders = allSurroundingTriangles(i,:); 
       % alle Ecken der umliegenden Dreiecke berechnen. Fallunterscheidung
       % ist notwendig, da 0 als Index nicht zugelassen ist. 
       alledges = zeros(6,3); 
       for j=1:6
           if surrounders(j) == 0 
               alledges(j,:) = 0; 
           else
                alledges(j,:) = edges(surrounders(j),:); 
           end
       end
       % hier will ich wieder mit Sparsematrizen arbeiten. Also betrachte
       % ich nur die Diagonalen, auf denen Werte stehen. Die Herleitung
       % dieser Formeln findet sich in der zugehörigen Bachelorarbeit. 
       % fn0 wird benötigt, da alledges auch 0 werden kann. Null ist als
       % Index nicht zugelassen. 
       % diag ist die Hauptdiagonale. Sie beinhaltet die Werte für T_i T_i
       diag(i) =  1/60 * (fn0(f,alledges(1,1)) + fn0(f,alledges(1,2)) + 3*fn0(f,alledges(1,3)) + fn0(f,alledges(2,1)) + fn0(f,alledges(2,2)) + 3*fn0(f,alledges(2,3)) + fn0(f,alledges(3,1)) + 3 * fn0(f,alledges(3,2)) + fn0(f,alledges(3,3)) + fn0(f,alledges(4,1)) + 3 * fn0(f,alledges(4,2)) + fn0(f,alledges(4,3)) + 3 * fn0(f,alledges(5,1)) + fn0(f,alledges(5,2)) + fn0(f,alledges(5,3)) + 3 * fn0(f,alledges(6,1)) + fn0(f,alledges(6,2)) + fn0(f,alledges(6,3))); 
       % secdiag ist die erste Nebendiagonale. Sie gibt an, wenn die
       % Gitterpunkte i und j direkt nebeneinander liegen, also T_i T_i+1
       secdiag(i) =   1/120 * (fn0(f,alledges(3,1)) + 2 * fn0(f,alledges(3,2)) + 2*fn0(f,alledges(3,3)) + 2*fn0(f,alledges(6,1)) + 2 * fn0(f,alledges(6,2)) +fn0(f,alledges(6,3)) ) ;
       % thirddiag ist die nächste Nebendiagonale. Sie liegt weiter weg, da
       % sie angibt, wenn j direkt unter dem Gitterpunkt i ist, also
       % j=i+n+1
       thirddiag(i) = 1/120 * (2*fn0(f,alledges(4,1)) + 2 * fn0(f,alledges(4,2)) + fn0(f,alledges(4,3)) + fn0(f,alledges(5,1)) + 2 * fn0(f,alledges(5,2)) + 2*fn0(f,alledges(5,3)) )