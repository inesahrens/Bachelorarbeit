%OPTIMAZATION berechnet die Optimierung von \int_\Omega (v^2 + eps_1) 
%|\nabla u|^2 + eps_2 |\nabla v|^2 + 1/eps_3
% (1-v)^2 dx
%   Dabei wurde die Optimierung in zwei Teile geteilt: einmal die
%   Optimierung nach u und die Optimierung nach v. 
%
%   Optimierung nach u: 
%   Die Optimierung nach u kann direkt gel�st werden mittels 
%   \int_\Omega (v^2 + eps_1) \nabla u_1
%   \nabla T_i dx \forall i \in N, wobei N die Anzahl der St�tzstellen ist.
%
%   Optimierung nach v:
%   Hier wird eine Semiglatte Newton Methode verwendet, da keine direkte 
%   L�sung m�glich ist.
% 
%   Insgesamt sieht das Verfahren wie folgt aus:
%   F�r k=1,2,... 
%   1. Berechne u durch l�sen eines Gleichungssystems
%   2. L�se das gleichungssystem 
%         -(G_1(v^k,eta^k)) = (G_1v  G_1eta) (s_1^k)
%          (G_2(v^k,eta^k))   (G_2v  G_2eta) (s_2^k)
%      nach s^k
%   3. v^k+1   = s_1^k + v^k
%      eta^k+1 = s_2^k + eta^k
%   4. wiederhole

%   Dabei ist G gegeben durch
%   G(v^k,eta^k)_1 = \int_\Omega T_i v |\nabla u|^2 + eps_2 \nabla v \nabla
%   T_i - 2/ eps_2 (1-v)T_i + eta T_i dx \forall i 
%   G(v^k,eta^k)_2 = \int_\Omega (eta - max{0,eta + c(v-v0)} 
%                       - min{0,eta+cv}) T_i dx \forall i 
%   wobei v_0 eine geg. Schranke ist und c ein Konstante.
%   Die genauen Ableitungen stehen an der Stelle, an der wir sie sp�ter
%   brauchen. 


% Deklaration der Variablen

% n+1 ist die Anzahl der St�tzstellen nach rechts
% m+1 ist die Anzahl der St�tzstellen nach unten
n=100;
m=100;
nu = .1; 
% Deklaration aller epsilons
epsilon = [0.01;nu* 5/(n+1) ;(5/(n+1) )*(1/nu)]; 

% setzten von v f�r Anfangsdaten. Hier liegt ein Riss in der Mitte des 
% Gebietes vor.  
v = ones((m+1)*(n+1),1); 
for i=0:20
    v(i*(n+1)+ n/2) = 0;
    v(i*(n+1)+ n/2+1) = 0;  
end

% eta ist der Lagrangemultiplikator. Er muss auch vorab gesetzt werden. 
eta = -.34*ones((n+1)*(m+1),1);  

%v0 ist die Schranke. Im Problem gilt, dass 0<v<v0
v0 =.9*ones((n+1)*(m+1),1); 

% die Konstante kommt durch den Lagrangemultiplikator rein. 
const = 1; 

%Anzahl der Newtonschritte die durchgef�hrt werden sollen
k=200;
%f�r plots
l=1;
% u0 sind die Randdaten von u. Hier wird u am rechten und am linken Rand
% eingespannt, sodass an einem Rand u 1 und und am anderen 2. 
u0=zeros((n+1)*(m+1), 1);

for i=0:m
     u0(i*(n+1)+1) = 1; 
     u0(i*(n+1)+n+1) = 2;
end

% Da alles mit dreieckig linearen Lagrangeelementen implementiert ist, muss
% das Gebiet zun�chst Triangulisiert werden. 
% edges ist eine n*m*2 x 3 Matrix und berechnet alle Ecken der Dreiecke
% allSurroundingTriangles ist eine (m+1)*(n+1) x 6 Matrix und beinhaltet in
% der i ten Zeile alle indizes der Dreiecke, die um den Punkt i liegen. 
[edges, allSurroundingTriangles] = triangulasation(m,n); 

for i=1:k

    % u1 wird berechnet. 
    u1 = calculate_u(edges, allSurroundingTriangles, u0, epsilon(1), v,...
        m, n);  
    u = [u1,u1] ; 

    % als n�chstes soll das Gleichungssystem gel�st werden, Dazu muss G
    % berechnte werden. Alle Matrizen, die gebraucht werden, werden hier
    % berechnet. 
    % numerisch dargestellt ist G1 nichts anderes als 
    % (A+2eps_2B+2/eps_3D)v-2/eps_3*D*ones() +D*eta
    % mit A = \int\Omega |\nabla u|^2 T_i T_j dx
    %     B = \int\Omega \nabla T_i \nabla T_j dx
    %     D = \int\Omega T_i T_j dx
    [A, B, D] = Matrixes_for_G(u, edges, allSurroundingTriangles, m, n ); 
    D=1/(n*m)*D; 
    % Berechnung von G1
    G1 = (2*A+2*epsilon(2)*B+ (2/epsilon(3))*D)*v - ...
        (2/epsilon(3))*D*ones((n+1)*(m+1),1) + D*eta;
    
    % Berechnung von G1v. Die Ableitung ist recht einfach zu sehen
    G1v = m*A + 2*epsilon(2)*B + (2/epsilon(3))*D; 

    % Berechnung von G1eta. Die Ableitung ist recht einfach zu sehen
    G1eta = D; 
    
    % Berechnung von G2 und den Ableitungen
    [G2, G2v, G2eta] = G2_and_deratives(v, eta, const, v0, D, edges,...
        allSurroundingTriangles, m, n); 
    
    % Die Ableitungen m�ssen noch zu einer gro�en Matrix, wie oben
    % beschrieben zusammengef�gt werden. Hier bieten sich sparse Matrizen
    % an, da es viele Nullen als Eintr�ge gibt. 
    gradG = sparse ((m+1)*(n+1)*2 ,(m+1)*(n+1)*2 ) ; 
    
    % G1v steht oben links
    gradG(1:(m+1)*(n+1),1:(m+1)*(n+1)) = G1v; 
    % G1eta steht oben rechts
    gradG(1:(m+1)*(n+1),(m+1)*(n+1)+1:(m+1)*(n+1)*2) = G1eta; 
    % G2v steht unten links
    gradG((m+1)*(n+1)+1:(m+1)*(n+1)*2,1:(m+1)*(n+1)) = G2v;     
    % G2eta steht unten rechts
    gradG((m+1)*(n+1)+1:(m+1)*(n+1)*2,(m+1)*(n+1)+1:(m+1)*(n+1)*2) = G2eta; 
    
    % G1 und G2 m�ssen als Vektor dargestllt werden. 
    G = sparse ((m+1)*(n+1)*2,1) ;
    % G1 steht oben
    G(1:(m+1)*(n+1),1) = - G1;
    % G2 steht unten
    G((m+1)*(n+1)+1:(m+1)*(n+1)*2,1) = - G2; 
    
    % Berechnet die Schrittweite s, indem das Gleichungssystem G = G's
    % gel�st wird. 
    s = gradG\G ;
    % berechnet das neue v
    v = s(1:(m+1)*(n+1)) + v;
    % berechnet das neue eta
    eta = s((m+1)*(n+1)+ 1:(m+1)*(n+1)*2) + eta;
    
    if i==1 || i==2 || i==5 || i==10 || i==20 || i==50 || i==100 || i==200
        
        resultu = reshape(u(:,1),n+1,m+1); 
        resultv = reshape(v,n+1,m+1); 
        resulteta = reshape(eta,n+1,m+1);
        [X,Y] = meshgrid(0:1/n:1);
        figure(3)
            subplot(2,4,l);
            mesh(X,Y,resultu);
        figure(4)  
            subplot(2,4,l);
            mesh(X,Y,resultv);
        l=l+1    
    end
end
