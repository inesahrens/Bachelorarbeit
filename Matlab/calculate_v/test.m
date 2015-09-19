%% Test Berechnung von v
% Testet alle Funktionen, die zur Berechnung von v notwendig sind. 

function result = test

% ResultGradUTotal = test_grad_u_total; 
% ResultGradUSurrounding = test_grad_u_surrounding; 
% ResultIntegrateGtigtj = test_integrate_gti_gtj; 
% ResultIntegratetitj = test_integrate_ti_tj; 
% ResultIntegrateUtitj = test_integrate_u_ti_tj; 
% ResultG2 = test_G_2; 
% ResultG2v = test_G_2_v; 
% ResultG2eta = test_G_2_eta; 
% ResultIntegrateEtaTiTriangles = test_integrate_eta_ti_triangles; 
% ResultIntegrateEtaTiTotal = test_integrate_eta_ti_total; 
% ResultIntegrateTi = test_integrate_ti
resultG2Deratives = test_g2_and_deratives

%% Section 0 Testet  G2_and_deratives(v, eta, const, v0, D, edges, allsurroundingTriangles, m, n)
% bla
function results = test_g2_and_deratives

testNum = 1; 

while true
    [v, eta, const, v0, D, edges, allSurroundingTriangles, m, n, correctG2, correctG2v, correctG2eta] = testCaseG2Deratives(testNum);
    
    if isnan(correctG2) % No test to run
        break
    end
    [codeG2, codeG2v, codeG2eta] = G2_and_deratives(v, eta, const, v0, D, edges, allSurroundingTriangles, m, n);
    
    results(testNum,1) = compareResults(codeG2, correctG2) ;  
    results(testNum,2) = compareResults(codeG2v, correctG2v) ; 
    results(testNum,3) = compareResults(codeG2eta, correctG2eta) ; 
    testNum = testNum + 1; 
end


function [v, eta, const, v0, D, edges, allSurroundingTriangles, m, n, correctG2, correctG2v, correctG2eta] = testCaseG2Deratives(testNum)
close all
switch testNum
    case 1
        %trivialer Fall
        m = 2; 
        n = 2; 
        [edges, allSurroundingTriangles] = triangulasation(m,n); 
        v = zeros((m+1)*(n+1),1); 
        eta = zeros((m+1)*(n+1),1);  
        const = 1; 
        v0 = ones((m+1)*(n+1),1); 
        u = zeros((n+1)*(m+1),2);
        [A, B, c, D] = Matrixes_for_G(u, edges, allSurroundingTriangles, m, n ); 
        correctG2 = zeros((n+1)*(m+1),1);
        correctG2v = zeros((n+1)*(m+1));
        correctG2eta = zeros((n+1)*(m+1));
    case 2
        % 1. Fall konstant
        m = 2; 
        n = 2; 
        [edges, allSurroundingTriangles] = triangulasation(m,n); 
        v = .5*ones((m+1)*(n+1),1); 
        eta = ones((m+1)*(n+1),1);  
        const = 1; 
        v0 = ones((m+1)*(n+1),1); 
        u = zeros((n+1)*(m+1),2);
        [A, B, c, D] = Matrixes_for_G(u, edges, allSurroundingTriangles, m, n ); 
        correctG2 = .5 * [1/3; 1/2;1/6; 1/2; 1; 1/2; 1/6; 1/2; 1/3];
        correctG2v = - D;
        correctG2eta = zeros((n+1)*(m+1));
    case 3
        % 2. Fall konstant
        m = 2; 
        n = 2; 
        [edges, allSurroundingTriangles] = triangulasation(m,n); 
        v = .5*ones((m+1)*(n+1),1); 
        eta = zeros((m+1)*(n+1),1);   
        const = 1; 
        v0 = ones((m+1)*(n+1),1); 
        u = zeros((n+1)*(m+1),2);
        [A, B, c, D] = Matrixes_for_G(u, edges, allSurroundingTriangles, m, n ); 
        correctG2 = zeros((n+1)*(m+1),1);
        correctG2v = zeros((n+1)*(m+1));
        correctG2eta = D;
    case 4
        % 3. Fall konstant
        m = 2; 
        n = 2; 
        [edges, allSurroundingTriangles] = triangulasation(m,n); 
        v = .5*ones((m+1)*(n+1),1); 
        eta = .5*ones((m+1)*(n+1),1);   
        const = 1; 
        v0 = ones((m+1)*(n+1),1); 
        u = zeros((n+1)*(m+1),2);
        [A, B, c, D] = Matrixes_for_G(u, edges, allSurroundingTriangles, m, n ); 
        correctG2 = .5 * [1/3; 1/2;1/6; 1/2; 1; 1/2; 1/6; 1/2; 1/3];
        correctG2v = zeros((n+1)*(m+1)); %-const/2 * D
        correctG2eta = zeros((n+1)*(m+1));   %.5 * D     
    case 5
        % übergang zwischen den Fällen
        % G2 korrekt, nur darstllungsfehler
        m = 2; 
        n = 4; 
        [edges, allSurroundingTriangles] = triangulasation(m,n); 
        v = .5*ones((m+1)*(n+1),1); 
        eta = [1;0.5;0;-.5;-1;1;0.5;0;-.5;-1;1;0.5;0;-.5;-1];   
        const = 1; 
        v0 = ones((m+1)*(n+1),1); 
        u = zeros((n+1)*(m+1),2);
        [A, B, c, D] = Matrixes_for_G(u, edges, allSurroundingTriangles, m, n ); 
        correctG2 = [1/6; 3/16; -1/24; -11/48; -1/12;  1/4; 5/12; 0; -5/12; -1/4;  1/12; 11/48; 1/24; -3/16; -1/6]; 
        
        diagG2v = - [7/60;1/60;0;1/20;1/15; 11/60;1/15;0;1/15;11/60; 1/15;1/20;0;1/60;7/60]; 
        secdiagG2v = - [1/60;0;0;1/40;0; 1/24;0;0;1/24;0; 1/40;0;0;1/60;0]; 
        thirddiagG2v = - [1/40;1/60;0;1/60;1/40; 1/40;1/60;0;1/60;1/40; 0;0;0;0;0]; 
        forthdiagG2v = - [1/24;0;0;1/24;0; 1/24;0;0;1/24;0; 0;0;0;0;0]; 
        secdiagG2v2 =  [0; secdiagG2v(1:(n+1)*(m+1)-1)]; 
        % berechnet die n+1 te Nebendiagonale von u
        thirddiagG2v2 = [zeros(n+1,1) ; thirddiagG2v(1:(n+1)*(m+1)-(n+1))] ;
        % berechnet die n+2 te Nebendiagonale von u
        forthdiagG2v2 = [zeros(n+2,1) ; forthdiagG2v(1:(n+1)*(m+1)-(n+2))] ;
        correctG2v = spdiags([forthdiagG2v thirddiagG2v secdiagG2v diagG2v secdiagG2v2 thirddiagG2v2 forthdiagG2v2],[-n-2 -n-1 -1 0 1 n+1 n+2],(n+1)*(m+1),(n+1)*(m+1) ); 
              
        diagG2eta = [0; 1/20; 11/60; 1/60; 0;  0; 1/15; 11/30; 1/15; 0;  0; 1/60; 11/60; 1/20; 0]; 
        secdiagG2eta = [0;1/40;1/60;0;0; 0;1/24;1/24;0;0; 0;1/60;1/40;0;0]; 
        thirddiagG2eta = [0;1/60;1/20;1/60;0; 0;1/60;1/20;1/60;0; 0;0;0;0;0 ]; 
        forthdiagG2eta =  [0;1/24;1/24;0;0; 0;1/24;1/24;0;0; 0;0;0;0;0] ; 
        secdiagG2eta2 =  [0; secdiagG2eta(1:(n+1)*(m+1)-1)]; 
        % berechnet die n+1 te Nebendiagonale von u
        thirddiagG2eta2 = [zeros(n+1,1) ; thirddiagG2eta(1:(n+1)*(m+1)-(n+1))] ;
        % berechnet die n+2 te Nebendiagonale von u
        forthdiagG2eta2 = [zeros(n+2,1) ; forthdiagG2eta(1:(n+1)*(m+1)-(n+2))] ;
        correctG2eta = spdiags([forthdiagG2eta thirddiagG2eta secdiagG2eta diagG2eta secdiagG2eta2 thirddiagG2eta2 forthdiagG2eta2],[-n-2 -n-1 -1 0 1 n+1 n+2],(n+1)*(m+1),(n+1)*(m+1) ); 
        
    otherwise 
        m = nan; 
        n = nan;  
        edges = nan; 
        allSurroundingTriangles = nan; 
        v = nan; 
        eta = nan; 
        const = nan;
        v0 = nan; 
        D = nan; 
        correctG2 = nan;
        correctG2v = nan; 
        correctG2eta = nan; 
end












%% Section 1 Testet  grad_u_total(u, edges, m, n )
% grad_u_total berechnet |\nabla u|^2 für alle Dreiecke der Triangulierung. 
% Dabei ist  |\nabla u|^2 = (u31-u21)^2 + (u11-u21)^2 + (u32-u22)^2 + (u12-u22)^2 
function results = test_grad_u_total

testNum = 1; 

while true
    [u, edges, m, n, correctData] = testCaseGradU(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = grad_u_total(u, edges, m, n ); 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end


function [u, edges, m, n, correctData] = testCaseGradU(testNum)
close all
switch testNum
    case 1
        m = 2; 
        n = 2; 
        edges = edges_of_triangles(m,n); 
        u = ones((n+1)*(m+1),2); 
        correctData = zeros(2*n*m,1); 
    case 2
        m = 2; 
        n = 2;  
        edges = edges_of_triangles(m,n); 
        u(:,1)=[1;2;3;2;3;4;3;4;5];
        u(:,2)=[1;2;3;2;3;4;3;4;5];
        correctData = 4* ones(2*n*m,1) ; 
    case 3
        m = 4; 
        n = 4;  
        edges = edges_of_triangles(m,n); 
        u(:,1)=[1;2;3;4;5;2;3;4;5;6;3;4;5;6;7;4;5;6;7;8;5;6;7;8;9];
        u(:,2)=[1;2;3;4;5;2;3;4;5;6;3;4;5;6;7;4;5;6;7;8;5;6;7;8;9];
        correctData = 4* ones(2*n*m,1) ;   
    case 4
        m = 2; 
        n = 2;  
        edges = edges_of_triangles(m,n); 
        u(:,1)=[1;2;1;1;2;1;1;2;1];
        u(:,2)=[1;2;1;1;2;1;1;2;1];
        correctData =2*ones(2*n*m,1) ;               
    otherwise 
        m = nan; 
        n = nan;  
        edges = nan; 
        u = nan; 
        correctData = nan; 
end

%% Section 2 testet grad_u_surrounding
% grad_u_surrounding berechnet |nable u|^2 von allen umliegenden Dreiecken

function results = test_grad_u_surrounding

testNum = 1; 

while true
    [surroundings, gradUTotal, correctData] = testCaseGradUSurrounding(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = grad_u_surroundings(surroundings, gradUTotal); 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end

function [surroundings, gradUTotal, correctData] = testCaseGradUSurrounding(testNum)
close all
switch testNum
    case 1
        surroundings = [0 0 0 0 1 2] ; 
        gradUTotal = 2*ones(1,2);  
        correctData = [0 0 0 0 2 2]; 
    case 2
        surroundings = [1 2 3 6 7 8] ; 
        gradUTotal = 2*ones(1,9);    
        correctData = 2*ones(1,6) ;         
    otherwise 
        surroundings = nan ;
        gradUTotal = nan ;  
        correctData = nan;    
end


%% Section 3 testet  integrate_gti_gtj( allSurroundingTriangles, m, n)
% integrate_gti_gtj berechnet integral nabla Ti nabla Tj dx

function results = test_integrate_gti_gtj

testNum = 1; 

while true
    [allSurroundingTriangles, m, n, correctData] = testCaseIntegrategTigTj(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = integrate_gti_gtj(allSurroundingTriangles, m, n); 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end

function [allSurroundingTriangles, m, n, correctData] = testCaseIntegrategTigTj(testNum)
close all
switch testNum
    case 1
        m = 2;
        n = 2;
        allSurroundingTriangles = surrounding_triangles(m,n); 
        correctData =0.5*[2,-1,0,-1,0,0,0,0,0;-1,4,-1,0,-2,0,0,0,0;0,-1,2,0,0,-1,0,0,0;-1,0,0,4,-2,0,-1,0,0;0,-2,0,-2,8,-2,0,-2,0;0,0,-1,0,-2,4,0,0,-1;0,0,0,-1,0,0,2,-1,0;0,0,0,0,-2,0,-1,4,-1;0,0,0,0,0,-1,0,-1,2]  ; 
        %     case 2
%         m = 
%         n = 
%         allSurroundingTriangles = 
%         correctData = ;        
    otherwise 
        m = nan; 
        n = nan; 
        allSurroundingTriangles = nan; 
        correctData = nan;    
end
 
%% Section 4 testet  integrate_ti_tj( allSurroundingTriangles, m, n)
% integrate_ti_tj berechnet integral Ti Tj dx

function results = test_integrate_ti_tj

testNum = 1; 

while true
    [allSurroundingTriangles, m, n, correctData] = testCaseIntegrateTiTj(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData =  integrate_ti_tj(allSurroundingTriangles, m, n); 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end

function [allSurroundingTriangles, m, n, correctData] = testCaseIntegrateTiTj(testNum)
close all
switch testNum
    case 1
        m = 2; 
        n = 2; 
        allSurroundingTriangles = surrounding_triangles(m,n); 
        correctData = [1/6,1/24,0,1/24,1/12,0,0,0,0;1/24,1/4,1/24,0,1/12,1/12,0,0,0;0,1/24,1/12,0,0,1/24,0,0,0;1/24,0,0,1/4,1/12,0,1/24,1/12,0;1/12,1/12,0,1/12,1/2,1/12,0,1/12,1/12;0,1/12,1/24,0,1/12,1/4,0,0,1/24;0,0,0,1/24,0,0,1/12,1/24,0;0,0,0,1/12,1/12,0,1/24,1/4,1/24;0,0,0,0,1/12,1/24,0,1/24,1/6]; 
%     case 2
%         m = ; 
%         n = ; 
%         allSurroundingTriangles = ; 
%         correctData = ;        
    otherwise 
        m = nan; 
        n = nan; 
        allSurroundingTriangles = nan; 
        correctData = nan;    
end


%% Section 5 testet  integrate_ti( allSurroundingTriangles, m, n)
% integrate_ti berechnet integral Ti  dx

function results = test_integrate_ti

testNum = 1; 

while true
    [allSurroundingTriangles, m, n, correctData] = testCaseIntegrateTi(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData =  integrate_ti(allSurroundingTriangles, m, n);  
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end

function [allSurroundingTriangles, m, n, correctData] = testCaseIntegrateTi(testNum)
close all
switch testNum
    case 1
        m = 2; 
        n = 2; 
        allSurroundingTriangles = surrounding_triangles(m,n); 
        correctData = [1/3;1/2;1/6;1/2;1;1/2;1/6;1/2;1/3];         
    otherwise 
        m = nan; 
        n = nan; 
        allSurroundingTriangles = nan; 
        correctData = nan;    
end




%% Section 6 testet  integrate_u_ti_tj(gradUTotal, allSurroundingTriangles, m, n)
% integrate_u_ti_tj berechnet integral |nabla u|^2 Ti Tj dx

function results = test_integrate_u_ti_tj

testNum = 1; 

while true
    [gradUTotal, allSurroundingTriangles, m, n, correctData] = testCaseIntegrateuTiTj(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = integrate_u_ti_tj(gradUTotal, allSurroundingTriangles, m, n); 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end

function [gradUTotal, allSurroundingTriangles, m, n, correctData] = testCaseIntegrateuTiTj(testNum)
close all
switch testNum
    case 1
        gradUTotal = 2*ones(9,1); 
        m = 2; 
        n = 2; 
        allSurroundingTriangles = surrounding_triangles(m,n); 
        correctData = 2*[1/6,1/24,0,1/24,1/12,0,0,0,0;1/24,1/4,1/24,0,1/12,1/12,0,0,0;0,1/24,1/12,0,0,1/24,0,0,0;1/24,0,0,1/4,1/12,0,1/24,1/12,0;1/12,1/12,0,1/12,1/2,1/12,0,1/12,1/12;0,1/12,1/24,0,1/12,1/4,0,0,1/24;0,0,0,1/24,0,0,1/12,1/24,0;0,0,0,1/12,1/12,0,1/24,1/4,1/24;0,0,0,0,1/12,1/24,0,1/24,1/6]; 
%     case 2
%         gradUTotal = ; 
%         m = ; 
%         n = ; 
%         allSurroundingTriangles = ; 
%         correctData = ;        
    otherwise 
        gradUTotal = nan; 
        m = nan; 
        n = nan; 
        allSurroundingTriangles = nan; 
        correctData = nan;    
end


%% Section 7 testet  G_2_eta(v, eta, c, v0, m, n )
%G2 berechnet G2v(v,eta)= -c         falls -c(v-v0) < eta oder eta < -cv 
%                         0          falls -cv < eta <  -c(v-v0) 
%                         in [-c,0]  falls -c(v-v0) = eta oder eta = -cv 

function results = test_G_2_eta

testNum = 1; 

while true
    [v, eta, c, v0, m, n, correctData] = testCaseG2eta(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = G_2_eta(v, eta, c, v0, m, n ); 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end

function [v, eta, c, v0, m, n, correctData] = testCaseG2eta(testNum)
close all
switch testNum
    case 1
        % 3. Fall im Programm
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        eta = zeros(4,1); 
        c = 1; 
        v0 = zeros(4,1); 
        correctData = 0.5*ones(4,1); 
    case 2
        %2. Fall im Programm
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        eta = .5*ones(4,1); 
        c = 1; 
        v0 = ones(4,1); 
        correctData = ones(4,1);    
    case 3
        %1. Fall im programm
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        eta = [-1,-1,2,2]; 
        c = 1; 
        v0 = ones(4,1); 
        correctData = zeros(4,1);   
    otherwise 
        m = nan; 
        n = nan; 
        v = nan;
        eta = nan; 
        c = nan; 
        v0 = nan; 
        correctData = nan; 
end



%% Section 8 testet  G_2_v(v, eta, c, v0, m, n )
%G2 berechnet G2v(v,eta)= 0          falls -c(v-v0) < eta oder eta < -cv 
%                         1          falls -cv < eta <  -c(v-v0) 
%                         in [0,1]   falls -c(v-v0) = eta oder eta = -cv 

function results = test_G_2_v

testNum = 1; 

while true
    [v, eta, c, v0, m, n, correctData] = testCaseG2v(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = G_2_v(v, eta, c, v0, m, n ); 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end

function [v, eta, c, v0, m, n, correctData] = testCaseG2v(testNum)
close all
switch testNum
    case 1
        % 3. Fall im Programm
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        eta = zeros(4,1); 
        c = 1; 
        v0 = zeros(4,1); 
        correctData = -c/2*ones(4,1); 
    case 2
        %2. Fall im Programm
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        eta = .5*ones(4,1); 
        c = 1; 
        v0 = ones(4,1); 
        correctData = zeros(4,1);    
    case 3
        %1. Fall im programm
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        eta = [-1,-1,2,2]; 
        c = 1; 
        v0 = ones(4,1); 
        correctData = -c*ones(4,1);   
    otherwise 
        m = nan; 
        n = nan; 
        v = nan;
        eta = nan; 
        c = nan; 
        v0 = nan; 
        correctData = nan; 
end


%% Section 9 testet  G_2(v, eta, c, v0, m, n )
%G2 berechnet G2(v,eta)= -c(v-v0) falls -c(v-v0) <= eta
%                        eta      falls -cv < eta <  -c(v-v0) 
%                        -cv      falls  eta <=  -cv 

function results = test_G_2

testNum = 1; 

while true
    [v, eta, c, v0, m, n, correctData] = testCaseG2(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = G_2(v, eta, c, v0, m, n ) ; 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end

function [v, eta, c, v0, m, n, correctData] = testCaseG2(testNum)
close all
switch testNum
    case 1
        % 1. Fall im Programm
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        eta = zeros(4,1); 
        c = 1; 
        v0 = zeros(4,1); 
        correctData = zeros(4,1); 
    case 2
        %2. Fall im Programm
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        eta = .5*ones(4,1); 
        c = 1; 
        v0 = ones(4,1); 
        correctData = eta;    
    case 3
        %3. Fall im programm
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        eta = -1*ones(4,1); 
        c = 1; 
        v0 = ones(4,1); 
        correctData = zeros(4,1);   
    otherwise 
        m = nan; 
        n = nan; 
        v = nan;
        eta = nan; 
        c = nan; 
        v0 = nan; 
        correctData = nan; 
end


%% Section 10 Testet integrate_eta_ti_triangles(eta, edges, m, n )
%integrate_eta_ti_triangles berechnet für jedes Dreieck int eta t_i für alle Ti
%die an dem Dreieck anliegen.  
function results = test_integrate_eta_ti_triangles

testNum = 1; 

while true
    [eta, edges, m, n, correctData] = testCaseIntegrateEtaTiTriangles(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = integrate_eta_ti_triangles(eta, edges, m, n) ; 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end


function [eta, edges, m, n, correctData] = testCaseIntegrateEtaTiTriangles(testNum)
close all
switch testNum
    case 1
        m = 1; 
        n = 1; 
        edges = edges_of_triangles(m,n); 
        eta = ones((m+1)*(n+1),1); 
        correctData = 4 * ones(2*m*2,3); 
    case 2
        m = 2; 
        n = 2; 
        edges = edges_of_triangles(m,n); 
        eta = [1,2,3,1,2,3,1,2,3]; 
        correctData = [5,5,6;6,7,7;9,9,10;10,11,11;5,5,6;6,7,7;9,9,10;10,11,11];                       
    otherwise 
        m = nan; 
        n = nan;  
        edges = nan; 
        eta = nan; 
        correctData = nan; 
end


%% Section 2 Testet integrate_eta_ti_total(intEtaTi, sur, m, n )
%integrate_eta_ti_triangles berechnet für jedes Dreieck int eta t_i für alle Ti
%die an dem Dreieck anliegen.  
function results = test_integrate_eta_ti_total

testNum = 1; 

while true
    [intEtaTi, sur, m, n, correctData] = testCaseIntegrateEtaTiTotal(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = integrate_eta_ti_total(intEtaTi, sur, m, n) ; 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end


function [intEtaTi, sur, m, n , correctData] = testCaseIntegrateEtaTiTotal(testNum)
close all
switch testNum
    case 1
        m = 1; 
        n = 1; 
        sur = surrounding_triangles(m,n); 
        intEtaTi = 4*ones((m+1)*(n+1),3); 
        correctData = [8;4;4;8]; 
    case 2
        m = 2; 
        n = 2; 
        sur = surrounding_triangles(m,n); 
        intEtaTi = 4* ones((m+1)*(n+1),3); 
        correctData = [8;12;4;12;24;12;4;12;8];                       
    otherwise 
        m = nan; 
        n = nan;  
        sur = nan; 
        intEtaTi = nan; 
        correctData = nan; 
end



function result = compareResults(codeData,correctData)

result = isequal(codeData,correctData); 



