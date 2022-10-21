function AssembleGlobalStiffMatrix(obj)
K11 = obj.K11;
K13 = obj.K13;
K22 = obj.K22;
K23 = obj.K23;
K33 = obj.K33;
%> For term $\frac{\partial}{\partial x}\left (\frac{\partial u}{\partial y}\right )$
mesh = obj.meshUnion;
InnerEdge = mesh.InnerEdge;
BottomEdge = mesh.BottomEdge;
SurfaceBoundaryEdge = mesh.SurfaceBoundaryEdge;

cell = obj.meshUnion.cell;
K = mesh.K;  Np = cell.Np;
obj.StiffMatrix = zeros(K*Np);
%For this test, the mesh is uniform
Dx = diag(mesh.rx(:,1))*mesh.cell.Dr +  diag(mesh.sx(:,1))*mesh.cell.Ds +  diag(mesh.tx(:,1))*mesh.cell.Dt;
Dy = diag(mesh.ry(:,1))*mesh.cell.Dr +  diag(mesh.sy(:,1))*mesh.cell.Ds +  diag(mesh.ty(:,1))*mesh.cell.Dt;
Dz = diag(mesh.rz(:,1))*mesh.cell.Dr +  diag(mesh.sz(:,1))*mesh.cell.Ds +  diag(mesh.tz(:,1))*mesh.cell.Dt;
ElementMassMatrix = diag(mesh.J(:,1))*cell.M;

% For term $$k_{11}\frac{\partial v}{\partial x}\frac{\partial p}{\partial x} + k_{13}\frac{\partial v}{\partial \sigma}\frac{\partial p}{\partial x} +
% k_{22}\frac{\partial v}{\partial y}\frac{\partial p}{\partial y} + k_{23}\frac{\partial v}{\partial \sigma}\frac{\partial p}{\partial y} + \
% k_{31}\frac{\partial v}{\partial x}\frac{\partial p}{\partial \sigma} + k_{32}\frac{\partial v}{\partial y}\frac{\partial p}{\partial \sigma}\
% + k_{33}\frac{\partial v}{\partial \sigma}\frac{\partial p}{\partial \sigma}$$
for i = 1:K
    obj.StiffMatrix((i-1)*Np+1:i*Np,(i-1)*Np+1:i*Np) = -1 * (((diag(K11(:,i)) * Dx)' * ElementMassMatrix * Dx + (diag(K13(:,i)) * Dz)' * ElementMassMatrix * Dx + ...
        (diag(K22(:,i)) * Dy)' * ElementMassMatrix * Dy + (diag(K23(:,i)) * Dz)' * ElementMassMatrix * Dy + ...
        (diag(K13(:,i)) * Dx)' * ElementMassMatrix * Dz + (diag(K23(:,i)) * Dy)' * ElementMassMatrix * Dz + ...
        (diag(K33(:,i)) * Dz)' * ElementMassMatrix * Dz));
end

fm = repmat( (obj.meshUnion.InnerEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.InnerEdge.FToE(1,:)))', 1, 1 );
Tau = 1000 * bsxfun(@times,  (fm(:) )',...
    ( max(obj.meshUnion(1).cell.Nz, obj.meshUnion(1).cell.N) + 1 )*( max(obj.meshUnion(1).cell.Nz, obj.meshUnion(1).cell.N) + ...
    3 )/3 * obj.meshUnion(1).cell.Nface/2);
for face = 1:InnerEdge.Ne
    OP11 = zeros(Np);
    OP12 = zeros(Np);
    OP22 = zeros(Np);
    OP21 = zeros(Np);
    ele = InnerEdge.FToE(2*(face-1)+1);
    adjacentEle = InnerEdge.FToE(2*(face-1)+2);
    eidM = InnerEdge.FToN1(:,face);
    eidP = InnerEdge.FToN2(:,face);
    nx = InnerEdge.nx(:,face);
    ny = InnerEdge.ny(:,face);
    Js = InnerEdge.Js(:,face);
    FacialMass2d = diag(Js)*InnerEdge.M;
    %> For term $$\left \{k_{11}\frac{\partial v}{\partial x} \right \}[p]_x+  \left \{k_{13}\frac{\partial v}{\partial \sigma} \right \}[p]_x
    %  + \left \{ k_{22}\frac{\partial v}{\partial y} \right \}[p]_y + \left \{ k_{23}\frac{\partial v}{\partial \sigma} \right \}[p]_y$$.
    %  Here, both $v$ and $p$ are local variable
    TempDx11 = diag(K11(:,ele))*Dx;
    TempDz13 = diag(K13(:,ele))*Dz;
    TempDy22 = diag(K22(:,ele))*Dy;
    TempDz23 = diag(K23(:,ele))*Dz;
    
    OP11(:, eidM) = OP11(:, eidM) + 0.5*TempDx11(eidM, :)'*diag(nx)*FacialMass2d + 0.5*TempDz13(eidM, :)'*diag(nx)*FacialMass2d + ...
        0.5*TempDy22(eidM,:)'*diag(ny)*FacialMass2d + 0.5*TempDz23(eidM, :)'*diag(ny)*FacialMass2d; % vertical direction not consider here
    
    %> For term $$\left \{k_{11}\frac{\partial v}{\partial x} \right \}[p]_x+  \left \{k_{13}\frac{\partial v}{\partial \sigma} \right \}[p]_x
    %  + \left \{ k_{22}\frac{\partial v}{\partial y} \right \}[p]_y + \left \{ k_{23}\frac{\partial v}{\partial \sigma} \right \}[p]_y$$
    % Here, p is local while v is adjacent
    TempDx11 = diag(K11(:,adjacentEle))*Dx;
    TempDz13 = diag(K13(:,adjacentEle))*Dz;
    TempDy22 = diag(K22(:,adjacentEle))*Dy;
    TempDz23 = diag(K23(:,adjacentEle))*Dz;
    OP12(:, eidM) = OP12(:, eidM) + 0.5*TempDx11(eidP, :)'*diag(nx)*FacialMass2d + 0.5*TempDz13(eidP, :)'*diag(nx)*FacialMass2d + ...
        0.5*TempDy22(eidP,:)'*diag(ny)*FacialMass2d + 0.5*TempDz23(eidP, :)'*diag(ny)*FacialMass2d; % vertical direction not consider here
    
    %> For term $$\left \{k_{11}\frac{\partial v}{\partial x} \right \}[p]_x+  \left \{k_{13}\frac{\partial v}{\partial \sigma} \right \}[p]_x
    %  + \left \{ k_{22}\frac{\partial v}{\partial y} \right \}[p]_y + \left \{ k_{23}\frac{\partial v}{\partial \sigma} \right \}[p]_y$$
    % Here, both p and v are adjacent
    OP22(:, eidP) = OP22(:, eidP) + 0.5*TempDx11(eidP, :)'*diag(-1*nx)*FacialMass2d + 0.5*TempDz13(eidP, :)'*diag(-1*nx)*FacialMass2d + ...
        0.5*TempDy22(eidP,:)'*diag(-1*ny)*FacialMass2d + 0.5*TempDz23(eidP, :)'*diag(-1*ny)*FacialMass2d;
    
    %> For term $$\left \{k_{11}\frac{\partial v}{\partial x} \right \}[p]_x+  \left \{k_{13}\frac{\partial v}{\partial \sigma} \right \}[p]_x
    %  + \left \{ k_{22}\frac{\partial v}{\partial y} \right \}[p]_y + \left \{ k_{23}\frac{\partial v}{\partial \sigma} \right \}[p]_y$$
    % Here, p is adjacent while v is local
    TempDx11 = diag(K11(:,ele))*Dx;
    TempDz13 = diag(K13(:,ele))*Dz;
    TempDy22 = diag(K22(:,ele))*Dy;
    TempDz23 = diag(K23(:,ele))*Dz;
    
    OP21(:, eidP) = OP21(:, eidP) + 0.5*TempDx11(eidM, :)'*diag(-1*nx)*FacialMass2d + 0.5*TempDz13(eidM, :)'*diag(-1*nx)*FacialMass2d + ...
        0.5*TempDy22(eidM,:)'*diag(-1*ny)*FacialMass2d + 0.5*TempDz23(eidM, :)'*diag(-1*ny)*FacialMass2d;
    
    %> For term $$\left \{k_{11}\frac{\partial p}{\partial x}\right \}[v]_x+   \left \{k_{13}\frac{\partial p}{\partial \sigma}\right \}[v]_x+
    %  \left \{ k_{22}\frac{\partial p}{\partial y}\right \}[v]_y + \left \{ k_{23}\frac{\partial p}{\partial \sigma}\right \}[v]_y
    % $$. Here, both p and v are local
    
    TempDx11 = diag(K11(:,ele))*Dx;
    TempDz13 = diag(K13(:,ele))*Dz;
    TempDy22 = diag(K22(:,ele))*Dy;
    TempDz23 = diag(K23(:,ele))*Dz;
    
    OP11(eidM, :) = OP11(eidM, :) + 0.5 * diag(nx) * FacialMass2d * TempDx11(eidM, :) + 0.5 * diag(nx) * FacialMass2d * TempDz13(eidM, :)+...
        0.5*diag(ny)*FacialMass2d*TempDy22(eidM,:) + 0.5*diag(ny)*FacialMass2d*TempDz23(eidM, :);
    
    %> For term $$\left \{k_{11}\frac{\partial p}{\partial x}\right \}[v]_x+   \left \{k_{13}\frac{\partial p}{\partial \sigma}\right \}[v]_x+
    %  \left \{ k_{22}\frac{\partial p}{\partial y}\right \}[v]_y + \left \{ k_{23}\frac{\partial p}{\partial \sigma}\right \}[v]_y
    % $$. Here, p is local while v is adjacent
    
    OP12(eidP,:) = OP12(eidP, :) + 0.5 * diag(-1*nx)* FacialMass2d * TempDx11(eidM, :) + 0.5 * diag(-1*nx)* FacialMass2d * TempDz13(eidM, :) + ...
        0.5*diag(-1*ny)*FacialMass2d*TempDy22(eidM,:) + 0.5 * diag(-1*ny)*FacialMass2d*TempDz23(eidM,:);
    
    %> For term $$\left \{k_{11}\frac{\partial p}{\partial x}\right \}[v]_x+   \left \{k_{13}\frac{\partial p}{\partial \sigma}\right \}[v]_x+
    %  \left \{ k_{22}\frac{\partial p}{\partial y}\right \}[v]_y + \left \{ k_{23}\frac{\partial p}{\partial \sigma}\right \}[v]_y
    % $$. Here, both p and v are adjacent
    TempDx11 = diag(K11(:,adjacentEle))*Dx;
    TempDz13 = diag(K13(:,adjacentEle))*Dz;
    TempDy22 = diag(K22(:,adjacentEle))*Dy;
    TempDz23 = diag(K23(:,adjacentEle))*Dz;
    
    OP22(eidP, :) = OP22(eidP, :) + 0.5 * diag(-1*nx) * FacialMass2d * TempDx11(eidP, :) + 0.5*diag(-1*nx)*FacialMass2d*TempDz13(eidP, :)+...
        0.5*diag(-1*ny)*FacialMass2d*TempDy22(eidP,:) + 0.5*diag(-1*ny)*FacialMass2d*TempDz23(eidP, :);
    
    %> For term $$\left \{k_{11}\frac{\partial p}{\partial x}\right \}[v]_x+   \left \{k_{13}\frac{\partial p}{\partial \sigma}\right \}[v]_x+
    %  \left \{ k_{22}\frac{\partial p}{\partial y}\right \}[v]_y + \left \{ k_{23}\frac{\partial p}{\partial \sigma}\right \}[v]_y
    % $$. Here, p is adjacent while v is local
    
    OP21(eidM, :) = OP21(eidM, :) + 0.5 * diag(nx) * FacialMass2d *TempDx11(eidP,:) + 0.5*diag(nx)*FacialMass2d*TempDz13(eidP,:) + ...
        0.5*diag(ny)*FacialMass2d*TempDy22(eidP,:) + 0.5*diag(ny)*FacialMass2d*TempDz23(eidP,:);
    
    %> For term $-\tau [q_h]_x[v]_x$  check to here
    OP11(eidM, eidM) = OP11(eidM, eidM) - Tau(face) * FacialMass2d;
    OP12(eidP, eidM) = OP12(eidP, eidM) - Tau(face) * (-1) * FacialMass2d;
    OP22(eidP, eidP) = OP22(eidP, eidP) - Tau(face) * FacialMass2d;
    OP21(eidM, eidP) = OP21(eidM, eidP) - Tau(face) * (-1) * FacialMass2d;
    
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
        OP11;
    obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) + ...
        OP12;
    obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = ...
        obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        OP22;
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        OP21;
end

fm = repmat( (obj.meshUnion.BottomEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.BottomEdge.FToE(1,:)))', 1, 1 );
Tau = 1000*25*bsxfun(@times,  (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    3 )/3 * obj.meshUnion(1).cell.Nface/2);

for face = 1:BottomEdge.Ne
    OP11 = zeros(Np);
    OP12 = zeros(Np);
    OP22 = zeros(Np);
    OP21 = zeros(Np);
    ele = BottomEdge.FToE(2*(face-1)+1);
    adjacentEle = BottomEdge.FToE(2*(face-1)+2);
    eidM = BottomEdge.FToN1(:,face);
    eidP = BottomEdge.FToN2(:,face);
    nz = BottomEdge.nz(:,face);
    Js = BottomEdge.Js(:,face);
    FacialMass2d = diag(Js)*BottomEdge.M;
    
    % For term $$\left \{k_{31}\frac{\partial v}{\partial x} \right \}[p]_{\sigma}+ \left \{ k_{32}\frac{\partial v}{\partial y} \right
    % \}[p]_{\sigma} + \left \{ k_{33}\frac{\partial v}{\partial \sigma} \right \}[p]_{\sigma}$$.
    % Here, both v and p are local
    
    TempDx31 = diag(K13(:,ele))*Dx;
    TempDy32 = diag(K23(:,ele))*Dy;
    TempDz33 = diag(K33(:,ele))*Dz;
    
    OP11(:, eidM) = OP11(:, eidM) + 0.5*TempDx31(eidM, :)'*diag(nz)*FacialMass2d + 0.5 * TempDy32(eidM, :)'*diag(nz)*FacialMass2d + ...
        0.5 * TempDz33(eidM,:)'*diag(nz)*FacialMass2d; % Horizontal direction not consider here
    
    % For term $$\left \{k_{31}\frac{\partial v}{\partial x} \right \}[p]_{\sigma}+ \left \{ k_{32}\frac{\partial v}{\partial y} \right
    % \}[p]_{\sigma} + \left \{ k_{33}\frac{\partial v}{\partial \sigma} \right \}[p]_{\sigma}$$.
    % Here, p is local while v is adjacent
    
    TempDx31 = diag(K13(:,adjacentEle))*Dx;
    TempDy32 = diag(K23(:,adjacentEle))*Dy;
    TempDz33 = diag(K33(:,adjacentEle))*Dz;
    
    OP12(:, eidM) = OP12(:, eidM) + 0.5*TempDx31(eidP,:)'*diag(nz)*FacialMass2d + 0.5 * TempDy32(eidP, :)'*diag(nz)*FacialMass2d + ...
        0.5 * TempDz33(eidP,:)'*diag(nz)*FacialMass2d;
    
    % For term $$\left \{k_{31}\frac{\partial v}{\partial x} \right \}[p]_{\sigma}+ \left \{ k_{32}\frac{\partial v}{\partial y} \right
    % \}[p]_{\sigma} + \left \{ k_{33}\frac{\partial v}{\partial \sigma} \right \}[p]_{\sigma}$$.
    % Here, both p and v are adjacent
    
    OP22(:, eidP) = OP22(:, eidP) + 0.5*TempDx31(eidP, :)'*diag(-1*nz)*FacialMass2d + 0.5 * TempDy32(eidP, :)'*diag(-1*nz)*FacialMass2d + ...
        0.5 * TempDz33(eidP,:)'*diag(-1*nz)*FacialMass2d;
    
    % For term $$\left \{k_{31}\frac{\partial v}{\partial x} \right \}[p]_{\sigma}+ \left \{ k_{32}\frac{\partial v}{\partial y} \right
    % \}[p]_{\sigma} + \left \{ k_{33}\frac{\partial v}{\partial \sigma} \right \}[p]_{\sigma}$$.
    % Here, here p is adjacent while v is local
    
    TempDx31 = diag(K13(:,ele))*Dx;
    TempDy32 = diag(K23(:,ele))*Dy;
    TempDz33 = diag(K33(:,ele))*Dz;
    
    OP21(:, eidP) = OP21(:, eidP) + 0.5*TempDx31(eidM, :)'*diag(-1*nz)*FacialMass2d + 0.5 * TempDy32(eidM, :)'*diag(-1*nz)*FacialMass2d + ...
        0.5 * TempDz33(eidM,:)'*diag(-1*nz)*FacialMass2d;
    
    %> For term $$\left \{k_{31}\frac{\partial p}{\partial x}\right \}[v]_{\sigma}+ \left \{ k_{32}\frac{\partial p}{\partial y}\right
    % \}[v]_{\sigma} + \left \{ k_{33}\frac{\partial p}{\partial \sigma}\right \}[v]_{\sigma}$$
    % Here, both p and v are local
    
    OP11(eidM, :) = OP11(eidM, :) + 0.5*diag(nz)*FacialMass2d*TempDx31(eidM, :) + 0.5*diag(nz)*FacialMass2d*TempDy32(eidM, :) + ...
        0.5 * diag(nz) * FacialMass2d*TempDz33(eidM, :);
    
    %> For term $$\left \{k_{31}\frac{\partial p}{\partial x}\right \}[v]_{\sigma}+ \left \{ k_{32}\frac{\partial p}{\partial y}\right
    % \}[v]_{\sigma} + \left \{ k_{33}\frac{\partial p}{\partial \sigma}\right \}[v]_{\sigma}$$
    % Here, p is local while v is adjacent
    
    OP12(eidP, :) = OP12(eidP, :) + 0.5*diag(-1*nz)*FacialMass2d*TempDx31(eidM, :) + 0.5*diag(-1*nz)*FacialMass2d*TempDy32(eidM, :) + ...
        0.5 * diag(-1*nz) * FacialMass2d*TempDz33(eidM, :);
    
    TempDx31 = diag(K13(:,adjacentEle))*Dx;
    TempDy32 = diag(K23(:,adjacentEle))*Dy;
    TempDz33 = diag(K33(:,adjacentEle))*Dz;
    
    %> For term $$\left \{k_{31}\frac{\partial p}{\partial x}\right \}[v]_{\sigma}+ \left \{ k_{32}\frac{\partial p}{\partial y}\right
    % \}[v]_{\sigma} + \left \{ k_{33}\frac{\partial p}{\partial \sigma}\right \}[v]_{\sigma}$$
    % Here, both p and v are adjacent
    
    OP22(eidP, :) = OP22(eidP, :) + 0.5*diag(-1*nz)*FacialMass2d*TempDx31(eidP, :) + 0.5*diag(-1*nz)*FacialMass2d*TempDy32(eidP, :) + ...
        0.5 * diag(-1*nz) * FacialMass2d*TempDz33(eidP, :);
    
    %> For term $$\left \{k_{31}\frac{\partial p}{\partial x}\right \}[v]_{\sigma}+ \left \{ k_{32}\frac{\partial p}{\partial y}\right
    % \}[v]_{\sigma} + \left \{ k_{33}\frac{\partial p}{\partial \sigma}\right \}[v]_{\sigma}$$
    % Here, p is adjacent while v is local
    
    OP21(eidM, :) = OP21(eidM, :) + 0.5*diag(nz)*FacialMass2d*TempDx31(eidP, :) + 0.5*diag(nz)*FacialMass2d*TempDy32(eidP, :) + ...
        0.5 * diag(nz) * FacialMass2d*TempDz33(eidP, :);
    
    %> For term $-\tau [q_h]_x[v]_x$  check to here
    OP11(eidM, eidM) = OP11(eidM, eidM) - Tau(face) * FacialMass2d;
    OP12(eidP, eidM) = OP12(eidP, eidM) - Tau(face) * (-1) * FacialMass2d;
    OP22(eidP, eidP) = OP22(eidP, eidP) - Tau(face) * FacialMass2d;
    OP21(eidM, eidP) = OP21(eidM, eidP) - Tau(face) * (-1) * FacialMass2d;
    
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
        OP11;
    obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) + ...
        OP12;
    obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = ...
        obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        OP22;
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        OP21;
end

fm = repmat( (obj.meshUnion.SurfaceBoundaryEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.SurfaceBoundaryEdge.FToE(1,:)))', 1, 1 );
Tau = 25 * 1000 * bsxfun(@times,  (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    3 )/3 * obj.meshUnion(1).cell.Nface/2);

if strcmp(obj.SurfaceBoundaryEdgeType, 'Dirichlet')
    for face = 1:SurfaceBoundaryEdge.Ne
        OP11 = zeros(Np);
        ele = SurfaceBoundaryEdge.FToE(2*(face-1)+1);
        eidM = SurfaceBoundaryEdge.FToN1(:,face);
        nz = SurfaceBoundaryEdge.nz(:,face);
        Js = SurfaceBoundaryEdge.Js(:,face);
        FacialMass2d = diag(Js)*SurfaceBoundaryEdge.M;
        
        
        % For term $$k_{31}\frac{\partial v}{\partial x} pn_{\sigma}+ k_{32}\frac{\partial v}{\partial y}pn_{\sigma} + k_{33}\frac{\partial v}{\partial \sigma} pn_{\sigma}$$
        % Here both p and v are local.
        TempDx31 = diag(K13(:,ele))*Dx;
        TempDy32 = diag(K23(:,ele))*Dy;
        TempDz33 = diag(K33(:,ele))*Dz;
        
        OP11(:, eidM) = OP11(:, eidM) + TempDx31(eidM, :)'*diag(nz)*FacialMass2d + TempDy32(eidM, :)'*diag(nz)*FacialMass2d + ...
            TempDz33(eidM,:)'*diag(nz)*FacialMass2d; % Horizontal direction not consider here
        
        % For term $$k_{31}\frac{\partial p}{\partial x} vn_{\sigma}+ k_{32}\frac{\partial p}{\partial y}vn_{\sigma} + k_{33}\frac{\partial p}{\partial \sigma} vn_{\sigma}$$
        % Here both p and v are local.
        
        OP11(eidM, :) = OP11(eidM, :) + diag(nz) * FacialMass2d * TempDx31(eidM, :) + diag(nz)*FacialMass2d*TempDy32(eidM, :)+...
            diag(nz)*FacialMass2d*TempDz33(eidM,:);
        %> For term $-\tau n_x^2q_h v$
        OP11(eidM,eidM) = OP11(eidM,eidM) - Tau(face)*FacialMass2d;
        obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
            OP11;
        
    end
elseif strcmp(obj.SurfaceBoundaryEdgeType, 'Newmann')

end
obj.StiffMatrix = sparse(obj.StiffMatrix);
end
