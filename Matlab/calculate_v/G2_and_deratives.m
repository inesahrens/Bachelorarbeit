function [ G2, G2v, G2eta ] = G2_and_deratives(v, eta, const, v0, D, edges, allSurroundingTriangles, m, n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    G2 = calculate_G2(v, eta, const, v0,D, m, n); 
    [f_v, f_eta] = calculate_f_v_eta(v, eta, const, v0, m, n); 
    
    G2v = calculate_G2div(f_v, edges, allSurroundingTriangles, m, n); 
    G2eta = calculate_G2div(f_eta, edges, allSurroundingTriangles, m, n); 
end
  
function G2 = calculate_G2(v, eta, const, v0, D, m, n)
%w berechnet w(v,eta)= -c(v-v0) falls -c(v-v0) <= eta
%                        eta      falls -cv < eta <  -c(v-v0) 
%                        -cv      falls  eta <=  -cv 
    w = zeros((n+1)*(m+1),1); 
    for i=1:(n+1)*(m+1)
        if ( -const*(v(i)-v0(i)) <= eta(i) ) 
            w(i) = -const*(v(i)-v0(i)) ; 
        elseif ( - const * v(i) < eta(i) && eta(i) <  -const*(v(i)-v0(i)) )
            w(i) = eta(i); 
        elseif ( eta(i) <= -const* v(i) )
            w(i) = - const * v(i) ; 
        end
    end
    G2 = D * w; 

end

function [f_v, f_eta] = calculate_f_v_eta(v, eta, const, v0, m, n)
%f_v berechnet f_v=      -c         falls -c(v-v0) < eta oder eta < -cv 
%                         0          falls -cv < eta <  -c(v-v0) 
%                         in [-c,0]  falls -c(v-v0) = eta oder eta = -cv 

%f_eta berechnet f_eta=   0          falls -c(v-v0) < eta oder eta < -cv 
%                         1          falls -cv < eta <  -c(v-v0) 
%                         in [0,1]   falls -c(v-v0) = eta oder eta = -cv 
    f_v = zeros((n+1)*(m+1),1); 
    f_eta = zeros((n+1)*(m+1),1); 
    
    for i=1:(n+1)*(m+1)
        for j=i:(n+1)*(m+1)
            if ( -const*(v(i)-v0(i)) < eta(i) || eta(i) < -const* v(i)) 
                f_v(i) = -const ; 
                f_eta(i) = 0 ; 
            elseif ( - const * v(i) < eta(i) && eta(i) <  -const*(v(i)-v0(i)) )
                f_v(i) = 0; 
                f_eta(j) = 1; 
            elseif ( -const*(v(i)-v0(i)) == eta(i) || eta(i) == -const* v(i) )
%                 f_v(i) = -const/2 ; 
%                 f_eta(i) =1/2; 
                  f_v(i) = 0; 
                  f_eta(i) = 0; 
            end
        end
    end  
end

function G2div = calculate_G2div(f, edges, allSurroundingTriangles, m, n)
    
    diag = zeros((n+1)*(m+1),1); 
    secdiag = zeros((n+1)*(m+1),1); 
    thirddiag = zeros((n+1)*(m+1),1);
    forthdiag = zeros((n+1)*(m+1),1); 
    for i=1:(n+1)*(m+1)
       surrounders = allSurroundingTriangles(i,:); 
       alledges = zeros(6,3); 
       for j=1:6
           if surrounders(j) == 0 
               alledges(j,:) = 0; 
           else
                alledges(j,:) = edges(surrounders(j),:); 
           end
       end
       diag(i) =  1/60 * (fn0(f,alledges(1,1)) + fn0(f,alledges(1,2)) + 3*fn0(f,alledges(1,3)) + fn0(f,alledges(2,1)) + fn0(f,alledges(2,2)) + 3*fn0(f,alledges(2,3)) + fn0(f,alledges(3,1)) + 3 * fn0(f,alledges(3,2)) + fn0(f,alledges(3,3)) + fn0(f,alledges(4,1)) + 3 * fn0(f,alledges(4,2)) + fn0(f,alledges(4,3)) + 3 * fn0(f,alledges(5,1)) + fn0(f,alledges(5,2)) + fn0(f,alledges(5,3)) + 3 * fn0(f,alledges(6,1)) + fn0(f,alledges(6,2)) + fn0(f,alledges(6,3))); 
       secdiag(i) =   1/120 * (fn0(f,alledges(3,1)) + 2 * fn0(f,alledges(3,2)) + 2*fn0(f,alledges(3,3)) + 2*fn0(f,alledges(6,1)) + 2 * fn0(f,alledges(6,2)) +fn0(f,alledges(6,3)) ) ;
       thirddiag(i) = 1/120 * (2*fn0(f,alledges(4,1)) + 2 * fn0(f,alledges(4,2)) + fn0(f,alledges(4,3)) + fn0(f,alledges(5,1)) + 2 * fn0(f,alledges(5,2)) + 2*fn0(f,alledges(5,3)) ) ;
       forthdiag(i) = 1/120 * (2 * fn0(f,alledges(5,1)) + fn0(f,alledges(5,2)) + 2 * fn0(f,alledges(5,3)) + 2 * fn0(f,alledges(6,1)) + fn0(f,alledges(6,2)) +2 * fn0(f,alledges(6,3)) ) ; 
    end
    % damit spDiags auf die richtigen Elemente zugreift, muss der erste
    % wert 0 gesetzt werden. 
    secdiag2 =  [0; secdiag(1:(n+1)*(m+1)-1)]; 
    % berechnet die n+1 te Nebendiagonale von u
    thirddiag2 = [zeros(n+1,1) ; thirddiag(1:(n+1)*(m+1)-(n+1))] ;
    % berechnet die n+2 te Nebendiagonale von u
    forthdiag2 = [zeros(n+2,1) ; forthdiag(1:(n+1)*(m+1)-(n+2))] ;
    G2div = spdiags([forthdiag thirddiag secdiag diag secdiag2 thirddiag2 forthdiag2],[-n-2 -n-1 -1 0 1 n+1 n+2],(n+1)*(m+1),(n+1)*(m+1) ); 
end

function result = fn0 (f,int)
    if int==0
        result = 0;
    else 
        result = f(int);
    end
end