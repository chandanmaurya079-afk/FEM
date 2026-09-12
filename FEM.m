clear; 
clc; 
close all;

% Data
ri     = 0.5;          
ro     = 1.0;        
k      = 10.0;         
h      = 25.0;         
Tinf   = 100.0;        
Tinner = 400.0;        
Q      = 0.0;          

nr     = 5;            
ntheta = 32;          

% Mesh generation 
rvec  = linspace(ri, ro, nr);
theta = linspace(0, 2*pi, ntheta+1); 
theta(end) = [];  % remove duplicate

nodes = zeros(nr*ntheta,2);
nid   = zeros(nr, ntheta);
id = 0;
for i = 1:nr
    for j = 1:ntheta
        id = id + 1;
        nodes(id,1) = rvec(i)*cos(theta(j));
        nodes(id,2) = rvec(i)*sin(theta(j));
        nid(i,j)    = id;
    end
end
nNodes = size(nodes,1);

% Connectivity
elems = [];
for i = 1:nr-1
    for j = 1:ntheta
        jp = j + 1; 
        if jp > ntheta, 
            jp = 1;
        end
        elems = [elems; nid(i,j) nid(i,jp) nid(i+1,jp) nid(i+1,j)];
    end
end
nElem = size(elems,1);

function [N, dN_dxi] = shapeQ4(xi, eta)
    
    % Shape functions (1x4 vector)
     N1 = 0.25*(1-xi)*(1-eta);
     N2 = 0.25*(1+xi)*(1-eta);
     N3 = 0.25*(1+xi)*(1+eta);
     N4 = 0.25*(1-xi)*(1+eta);
    
    N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta),(1+xi)*(1+eta), (1-xi)*(1+eta)];

    % Derivatives w.r.t. xi and eta 
    dN_dxi = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta); -(1-xi),  -(1+xi),  (1+xi),   (1-xi)];

end

%% Global matrices 
K = sparse(nNodes,nNodes);
F = zeros(nNodes,1);

% 2x2 Gauss
g = 1/sqrt(3);
gp = [-g -g; g -g; g g; -g g];
w  = [1 1 1 1];

% Element loop 
for e = 1:nElem
    conn = elems(e,:);
    Xe = nodes(conn,:);
    Ke = zeros(4,4); 
    fe = zeros(4,1);
    for kq = 1:4
        xi = gp(kq,1); 
        eta = gp(kq,2);
        [N, dN_dxi] = shapeQ4(xi, eta);
        J = dN_dxi * Xe;
        detJ = det(J);
        invJ = inv(J);
        dN_dx = invJ * dN_dxi;      % 2x4
        Ke = Ke + k * (dN_dx' * dN_dx) * detJ * w(kq);
        if Q ~= 0
            fe = fe + N' * Q * detJ * w(kq);
        end
    end
    K(conn,conn) = K(conn,conn) + Ke;
    F(conn)      = F(conn) + fe;
end

% Convection on outer boundary
gp1d = [-1/sqrt(3), 1/sqrt(3)];
w1d  = [1, 1];
for j = 1:ntheta
    jp = j + 1; if jp > ntheta, jp = 1; 
    end
    n1 = nid(nr,j); 
    n2 = nid(nr,jp);
    x1 = nodes(n1,1);
    y1 = nodes(n1,2);
    x2 = nodes(n2,1); 
    y2 = nodes(n2,2);
    L  = sqrt((x2-x1)^2 + (y2-y1)^2);
    Je = L/2;
    Ke_edge = zeros(2,2); 
    fe_edge = zeros(2,1);
    for gq = 1:2
        s = gp1d(gq);
        N1 = (1-s)/2;
        N2 = (1+s)/2;
        Nedge = [N1 N2];
        Ke_edge = Ke_edge + (Nedge'*h*Nedge) * Je * w1d(gq);
        fe_edge = fe_edge + (Nedge'*h*Tinf)  * Je * w1d(gq);
    end
    idx = [n1 n2];
    K(idx,idx) = K(idx,idx) + Ke_edge;
    F(idx)     = F(idx)     + fe_edge;
end

%Dirichlet inner boundary 
fixed = nid(1,:);
T = zeros(nNodes,1);
T(fixed) = Tinner;

F = F - K*T;
K(fixed,:) = 0; K(:,fixed)=0;
for ii = fixed
    K(ii,ii) = 1; 
    F(ii) = T(ii);
end

% Solve
T = K\F;

% Results 
fprintf('Min T = %.2f K,  Max T = %.2f K\n', min(T), max(T));

% Convert quads to triangles for plotting
tri = zeros(2*nElem,3);
for e = 1:nElem
    c = elems(e,:);
    tri(2*e-1,:) = c([1 2 3]);
    tri(2*e,:)   = c([1 3 4]);
end

figure('Color','w');
trisurf(tri, nodes(:,1), nodes(:,2), T, 'EdgeColor','none');
view(2); axis equal; colorbar;
title('Temperature Distribution (K)');
xlabel('x (m)'); ylabel('y (m)');
