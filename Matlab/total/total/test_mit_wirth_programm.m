
function test

  n = 11;
  n1=10;
  m=10; 
  
  x = linspace(0,1,n);
  [X,Y] = meshgrid(x,x);
  nodes = [X(:),Y(:)];
  elements = delaunay(nodes(:,1),nodes(:,2));

  u1 = zeros(size(X));
  u1(:,1:n/2) = 1;
  
  u = reshape(u1,n^2,1); 
  
  [edges, allSurroundingTriangles] = triangulasation(m,n1); 
  
Areal = weightedMassMatrix(nodes,elements,u1); 
[A, B, c, D] = Matrixes_for_G(u, edges, allSurroundingTriangles, m, n1 ); 



end



function M = weightedMassMatrix(nodes,elements,u)
% computes mass matrix M = (\int_\Omega |\nabla u|^2 T_i T_j dx)_{ij}
  n = length(nodes);
  rows = zeros(4*length(elements),1);
  cols = rows;
  valsM = rows;

  I = eye(2);
  D = [-1 -1;I]';

  counter = 1;
  for k = 1:size(elements,1)          % loop through all elements {T}
      nodeInds = elements(k, 1:3);    % get global indices
      x = nodes(nodeInds, 1:2);       % global nodes (as matrix rows)
      A = (D * x)';                   % Jacobian of coordinate transform from T_ref to T.
      AinvD = (D'/A)';                % gradient operator
      vol = abs(det(A)/2);            % triangle-volume
      gradU = AinvD*u(nodeInds);      % \nabla u

      % assemble matrix
      [c,r] = meshgrid(nodeInds,nodeInds);
      rows(counter:counter+8) = r(:);
      cols(counter:counter+8) = c(:);
      valsM(counter:counter+8) = [2,1,1;1,2,1;1,1,2] * (gradU'*gradU*vol/12);
      counter = counter + 9;
  end

  M = sparse(rows,cols,valsM,n,n);
end
