function result = test_u
% ruft alle Tests auf
    triangulasation = test_triangulasation %TODO: google correct spelling
    integral_of_v_total = test_integral_of_v_total
    integrate_u0 = test_integrate_u0
    integrate_u = test_integrate_u
end



%% Section 1 Testet  grad_u_total(u, edges, m, n )
% grad_u_total berechnet |\nabla u|^2 für alle Dreiecke der Triangulierung. 
% Dabei ist  |\nabla u|^2 = (u31-u21)^2 + (u11-u21)^2 + (u32-u22)^2 + (u12-u22)^2 
function results = test_triangulasation

testNum = 1; 

while true
    [m, n, correctE, correctEdge ] = testCaseTriangulasation(testNum);
    
    if isnan(correctE) % No test to run
        break
    end
    [codeEdge, codeE] = triangulasation(m,n); 
    results(testNum,1) = compareResults(codeE,correctE) ;  
    results(testNum,2) = compareResults(codeEdge,correctEdge) ; 
    testNum = testNum + 1; 
end
end

function [m, n, correctE, correctEdge] =  testCaseTriangulasation(testNum)
close all
switch testNum
    case 1
        m = 1; 
        n = 1; 
        correctE = [0,0,0,0,1,2; 0,0,0,2,0,0;0,0,1,0,0,0;1,2,0,0,0,0] ; 
        correctEdge = [1,3,4;1,2,4]; 
    case 2
        m = 2; 
        n = 2;  
        correctE = [0,0,0,0,1,2;0,0,0,2,3,4; 0,0,0,4,0,0; 0,0,1,0,5,6; 1,2,3,6,7,8; 3,4,0,8,0,0; 0,0,5,0,0,0; 5,6,7,0,0,0; 7,8,0,0,0,0]; 
        correctEdge = [1,4,5;1,2,5;2,5,6;2,3,6;4,7,8;4,5,8;5,8,9;5,6,9];          
    otherwise 
        m = nan; 
        n = nan;  
        correctE = nan; 
        correctEdge = nan; 
end
end


%% Section 2 testet integral_of_v_total(v, edges, epsilon, m, n)
%integral_of_v_total berechnet int v^2 + epsilon über alle dreiecke

function results = test_integral_of_v_total

testNum = 1; 

while true
    [v, edges, epsilon, m, n, correctData] = testCaseIntegralv(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = integral_of_v_total(v, edges, epsilon, m, n); 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end
end

function [v, edges, epsilon, m, n, correctData] = testCaseIntegralv(testNum)
close all
switch testNum
    case 1
        m = 1; 
        n = 1; 
        v = zeros(4,1);
        [edges, codeE] = triangulasation(m,n); 
        epsilon = 1; 
        correctData = epsilon * [1/2;1/2]; 
    case 2
        m = 1; 
        n = 1; 
        v = 2* ones(4,1);
        [edges, codeE] = triangulasation(m,n); 
        epsilon = 1; 
        correctData = 5 * [1/2;1/2];  
    case 3
        % rundungsfehler
        m = 1; 
        n = 1; 
        v = [1,1,2,2];
        [edges, codeE] = triangulasation(m,n); 
        epsilon = 1; 
        correctData = [17/12+1/2; 11/12+1/2]; 
    case 4
        m=2;
        n=2; 
        v=[1;0;1; 1;0;1; 1;0;1];
        [edges, codeE] = triangulasation(m,n); 
        epsilon = .5; 
        correctData =[1/4;1/12;1/12;1/4;1/4;1/12;1/12;1/4] + 1/4 *ones(8,1); 
    otherwise 
        m = nan; 
        n = nan; 
        v = nan;
        edges = nan; 
        epsilon = nan; 
        correctData = nan; 
end
end


%% Section 3 Testet integrate_u0(edges, allSurroundingTriangles, integral_v, u0,m,n )
%integrate_u0 berechnet \int_Omega (v^2 + eps) \nabla u0 \nabla T_j 

function results = test_integrate_u0

testNum = 1; 

while true
    [edges, allSurroundingTriangles, integral_v, u0, m, n, correctData] = testCaseintegrateU0(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData =  integrate_u0(edges, allSurroundingTriangles, integral_v, u0, m, n ); 
    results(testNum) = compareResults(codeData,correctData) ;  
    testNum = testNum + 1; 
end
end

function [edges, allSurroundingTriangles, integral_v, u0, m, n, correctData] = testCaseintegrateU0(testNum)
close all
switch testNum
    %alle tests (bis auf 1) laufen, wenn das ganze nicht vorher 0
    %gesetzt wird. 
    case 1
        m = 2;  
        n = 2;
        [edges, allSurroundingTriangles] = triangulasation(m, n); 
        epsilon = 1; 
        v = zeros(9,1);
        integral_v = integral_of_v_total(v, edges, epsilon, m, n); 
        u0 = zeros(9,1); 
        correctData = zeros(9,1); 
    case 2
        %geht schief, da Code davon ausgeht, dass anliegende neben neben
        %den rand 0 sind, ist aber hier nicht der fall
        m = 2; 
        n = 2;
        [edges, allSurroundingTriangles] = triangulasation(m,n); 
        epsilon = 1; 
        v = zeros(9,1); 
        integral_v = integral_of_v_total(v, edges, epsilon,m,n); 
        u0 = [1;0;1;1;0;1;1;0;1]; 
        correctData = .5*[1;-2;1;2;-4;2;1;-2;1];    
      case 3
        m = 1; 
        n = 3;
        [edges, allSurroundingTriangles] = triangulasation(m,n); 
        epsilon = 1; 
        v = zeros(8,1); 
        integral_v = integral_of_v_total(v, edges, epsilon,m,n); 
        u0 = [1;0;0;1;1;0;0;1]; 
        correctData = .5*[1;-1;-1;1;1;-1;-1;1];    
       case 4
        m = 2; 
        n = 3;
        [edges, allSurroundingTriangles] = triangulasation(m,n); 
        epsilon = 1; 
        v = zeros(12,1); 
        integral_v = integral_of_v_total(v, edges, epsilon,m,n); 
        u0 = [1;0;0;1;1;0;0;1;1;0;0;1]; 
        correctData = .5*[1;-1;-1;1; 2;-2;-2;2 ;1;-1;-1;1];    
    otherwise 
        m = nan; 
        n = nan;  
        edges = nan; 
        allSurroundingTriangles= nan; 
        integral_v = nan; 
        u0 = nan; 
        correctData = nan; 
end
end


%% Section 4 Testet  integrate_u(allSurroundingTriangles, integral_v, m, n)
% integrate_u berechnet das Integral über u, also \sum_i^{(n+1)(m+1)} \int_Omega (v^2+eps) \nabla T_i \nabla T_j

function results = test_integrate_u

testNum = 1; 

while true
    [allSurroundingTriangles, integral_v, m, n, correctData] = testCaseIntegrateU(testNum);
    
    if isnan(correctData) % No test to run
        break
    end
    codeData = integrate_u(allSurroundingTriangles, integral_v, m, n)  ; 
    b = full(codeData)
    results(testNum) = compareResults(codeData,correctData) ;   
    testNum = testNum + 1; 
end
end


function [allSurroundingTriangles, integral_v, m, n, correctData] =  testCaseIntegrateU(testNum)
close all
switch testNum
    case 1
        m = 1; 
        n = 1; 
        epsilon = 1;
        v = zeros(4,1); 
        [edges, allSurroundingTriangles] = triangulasation(m, n);  
        integral_v = integral_of_v_total(v, edges, epsilon, m, n) ; 
        correctData = .5*[2,-1,-1,0; -1,2,0,-1; -1,0,2,-1; 0,-1,-1,2]; 
    case 2
        % rundungsfehler
        m = 2; 
        n = 2; 
        epsilon = 0;
        v = [1,0,1,1,0,1,1,0,1]; 
        [edges, allSurroundingTriangles] = triangulasation(m, n);  
        integral_v = integral_of_v_total(v, edges, epsilon, m, n)  ; 
        correctData = [1/3,-1/12,0,-1/4,0,0,0,0,0; -1/12,1/2,-1/4,0,-1/6,0,0,0,0; 0,-1/4,1/2,0,0,-1/4,0,0,0; -1/4,0,0,5/6,-1/3,0,-1/4,0,0; 0,-1/6,0,-1/3,1,-1/3,0,-1/6,0;0,0,-1/4,0,-1/3,5/6,0,0,-1/4; 0,0,0,-1/4,0,0,1/2,-1/4,0; 0,0,0,0,-1/6,0,-1/4,1/2,-1/12; 0,0,0,0,0,-1/4,0,-1/12,1/3]         
    otherwise 
        m = nan; 
        n = nan;  
        allSurroundingTriangles = nan;
        integral_v = nan; 
        correctData = nan; 
end
end





function result = compareResults(codeData,correctData)

result = isequal(codeData,correctData); 
end