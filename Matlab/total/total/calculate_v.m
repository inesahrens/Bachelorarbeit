%OPTIMAZATION berechnet die Optimierung von \int_\Omega (v^2 + eps_1) |\nabla u|^2 + eps_2 |\nabla v|^2 + 1/eps_3
% (1-v)^2 dx
%   Dabei wurde die Optimierung in zwei Teile geteilt: einmal die
%   Optimierung nach u und die Optimierung nach v. 
%
%   Optimierung nach u: 
%   Die Optimierung nach u kann direkt gel�st werden mittels \int_\Omega (v^2 + eps_1) \nabla u_1
%   \nabla T_i dx \forall i \in N, wobei N die Anzahl der St�tzstellen ist.
%
%   Optimierung nach v:
%   Hier wird eine Semiglatte Newton Methode verwendet, da keine direkte L�sung m�glich ist.
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
%   G(v^k,eta^k)_2 = \int_\Omega (eta - max{0,eta + c(v-v0)} - min{0,eta+cv}) T_i dx \forall i 
%   wobei v_0 eine geg. Schranke ist und c ein Konstante.
%   Die genauen Ableitungen stehen an der Stelle, an der wir sie sp�ter
%   brauchen. 


% Deklaration der Variablen

% n+1 ist die Anzahl der St�tzstellen nach rechts
% m+1 ist die Anzahl der St�tzstellen nach unten
n=50;
m=50;
nu = .010; 
% Deklaration aller epsilons
epsilon = [0.01;nu* 20*1/n ;(20*1/n )/nu]; 

% setzten von v f�r Anfangsdaten. Hier liegt ein Riss in der Mitte des Gebietes vor.  
v = ones((m+1)*(n+1),1); 
for i=0:m
    v(i*(n+1)+n/2) = 0;
    v(i*(n+1)+n/2 + 1) = 0;
end


% eta ist der Lagrangemultiplikator. Er muss auch vorab gesetzt werden. 
eta = .5*ones((n+1)*(m+1),1);  

%v0 ist die Schranke. Im Problem gilt, dass 0<v<v0
v0 =2* ones((n+1)*(m+1),1); 

% die Konstante kommt durch den Lagrangemultiplikator rein. 
const = 1; 

%Anzahl der Newtonschritte die durchgef�hrt werden sollen
k=1; 

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
% allSurroundingTriangles ist eine (m+1)*(n+1) x 6 Matrix und beinhaltet in der
% i ten Zeile alle indizes der Dreiecke, die um den Punkt i liegen. 
[edges, allSurroundingTriangles] = triangulasation(m,n); 

u1 = calculate_u(edges, allSurroundingTriangles, u0, epsilon(1), v, m, n); 
u2 = calculate_u(edges, allSurroundingTriangles, u0, epsilon(1), v, m, n); 
u = [u1,u2] ; 


for i=1:k
%while (norm(v1-v)>.00002)

    %u1 und u2 bekommen eig unterschiedliche Randwerte? sonst sind sie
    %immer gleich? 
    % u1 und u2 werden berechnet. 
    %u1 = calculate_u(edges, allSurroundingTriangles, u0, epsilon(1), v, m, n); 
    %u2 = calculate_u(edges, allSurroundingTriangles, u0, epsilon(1), v, m, n); 
    
    
    %u1 = [1;1;2;2;2; 1;1;2;2;2; 1;1;2;2;2; 1;1;2;2;2; 1;1;2;2;2]; 
    %u2 = [1;1;2;2;2; 1;1;2;2;2; 1;1;2;2;2; 1;1;2;2;2; 1;1;2;2;2]; 
    %u = [u1,u2] ; 
    % als n�chstes soll das Gleichungssystem gel�st werden, Dazu muss G
    % berechnte werden. Alle Matrizen, die gebraucht werden, werden hier
    % berechnet. 
    % numerisch dargestellt ist G1 nichts anderes als 
    % (A+2eps_2B+2/eps_3D)v-2/eps_3c +e
    % mit A = \int\Omega |\nabla u|^2 T_i T_j dx
    %     B = \int\Omega \nabla T_i \nabla T_j dx
    %     c = \int\Omega T_i dx
    %     D = \int\Omega T_i T_j dx
    %     e = \int\Omega eta T_i dx
    [A, B, c, D] = Matrixes_for_G(u, edges, allSurroundingTriangles, m, n ); 
    e = integrate_eta_ti_total(edges, allSurroundingTriangles, eta, m, n ); 

    % Berechnung von G1
    G1 = (2*A+2*epsilon(2)*B+ (2/epsilon(3))*D)*v - (2/epsilon(3))*c + e; 
    
    % Berechnung von G1v. Die Ableitung ist recht einfach zu sehen
    G1v = 2*A + 2*epsilon(2)*B + (2/epsilon(3))*D; 

    % Berechnung von G1eta. Die Ableitung ist recht einfach zu sehen
    G1eta = D; 
    test = 4/epsilon(3) *G1v\(D*ones( (n+1)*(m+1),1))
    % Berechnung von G2 und den Ableitungen
    [G2, G2v, G2eta] = G2_and_deratives(v, eta, const, v0, D, edges, allSurroundingTriangles, m, n); 
    
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
    v = test; 
    % berechnet das neue eta
    eta = s((m+1)*(n+1)+ 1:(m+1)*(n+1)*2) + eta; 
end

% eigentlich sind u und v keine Vektoren, sondern Matrizen. 
resultu = reshape(u(:,1),n+1,m+1); 
resultv = reshape(v,n+1,m+1); 
resulteta = reshape(eta,n+1,m+1);


% Plotten des Ergebnisses 
% noch geht das nur f�r Gitter, die gleichhoch wie breit sind. 
[X,Y] = meshgrid(0:1/n:1);
figure 
mesh(X,Y,resultu)

figure
mesh(X,Y,resultv)

figure 
mesh(X,Y,resulteta)


