function Awirth = test_v_wirth

  % define mesh
  n = 5;
  x = linspace(0,1,n);
  [X,Y] = meshgrid(x,x);
  nodes = [X(:),Y(:)]; 
  elements = delaunay(nodes(:,1),nodes(:,2)); 
  [edges, allSurroundingTriangles] = triangulasation(n-1,n-1);
  elements = edges; 

  % define parameters
  nu = .1;
  epsilon = 5/n;
  
  % define u
  u = zeros(size(X));
  u(:,1:n/2) = 1;
  
  % compute v
  [M,L] = massStiffMatrix(nodes,elements);
  M2 = weightedMassMatrix(nodes,elements,u(:)); 
  A = M2 + nu * (epsilon*L+1/epsilon*M);
  Awirth = A; 
  rhs = M*ones(length(M),1)/epsilon;
  v = A\rhs;
  v = reshape(v,n,n);
  
  % visualisation
  subplot(1,2,1);
  surf(X,Y,u);
  subplot(1,2,2);
  surf(X,Y,v);

end

 
% nodes = n by 2 array whose rows are the mesh nodes
% elements = m by 3 array whose rows are the indices of the nodes of a triangle

function [M,L] = massStiffMatrix(nodes,elements)
% computes standard mass matrix M and stiffness matrix L using linear finite elements
  n = length(nodes);
  rows = zeros(4*length(elements),1);
  cols = rows;
  valsM = rows;
  valsL = rows;

  I = eye(2);
  D = [-1 -1;I]';

  counter = 1;
  for k = 1:size(elements,1)          % loop through all elements {T}
      nodeInds = elements(k, 1:3);    % get global indices
      x = nodes(nodeInds, 1:2);       % global nodes (as matrix rows)
      A = (D * x)';                   % Jacobian of coordinate transform from T_ref to T.
      AinvD = (D'/A)';                % gradient operator
      vol = abs(det(A)/2);            % triangle-volume

      % assemble matrices
      [c,r] = meshgrid(nodeInds,nodeInds);
      rows(counter:counter+8) = r(:);
      cols(counter:counter+8) = c(:);
      valsM(counter:counter+8) = [2,1,1;1,2,1;1,1,2] * (vol/12);
      valsL(counter:counter+8) = vol*AinvD'*AinvD;
      counter = counter + 9;
  end

  M = sparse(rows,cols,valsM,n,n);
  L = sparse(rows,cols,valsL,n,n);
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
