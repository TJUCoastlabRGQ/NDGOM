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

OP11 = -1 * (ElementMassMatrix\(K11 * Dx' * ElementMassMatrix * Dx + K13 * Dz' * ElementMassMatrix * Dx + ...
    K22 * Dy' * ElementMassMatrix * Dy + K23 * Dz' * ElementMassMatrix * Dy + ...
    K13 * Dx' * ElementMassMatrix * Dz + K23 * Dy' * ElementMassMatrix * Dz + ...
    K33 * Dz' * ElementMassMatrix * Dz));
for i = 1:K
    obj.StiffMatrix((i-1)*Np+1:i*Np,(i-1)*Np+1:i*Np) = OP11;
end
fm = repmat( (obj.meshUnion.InnerEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.InnerEdge.FToE(1,:)))', 1, 1 );
Tau = bsxfun(@times,  (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
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
    %> For term $[q_h]_y\left {\nabla_x v\right }$
    OP11(:, eidM) = OP11(:, eidM) + 0.5*K11*Dx(eidM, :)'*diag(nx)*FacialMass2d + 0.5*K13*Dz(eidM, :)'*diag(nx)*FacialMass2d + ...
                    0.5*K22*Dy(eidM,:)'*diag(ny)*FacialMass2d + 0.5*K23*Dz(eidM, :)'*diag(ny)*FacialMass2d; % vertical direction not consider here
    OP12(:, eidM) = OP12(:, eidM) + 0.5*K11*Dx(eidP, :)'*diag(nx)*FacialMass2d + 0.5*K13*Dz(eidP, :)'*diag(nx)*FacialMass2d + ...
                    0.5*K22*Dy(eidP,:)'*diag(ny)*FacialMass2d + 0.5*K23*Dz(eidP, :)'*diag(ny)*FacialMass2d; % vertical direction not consider here
    OP22(:, eidP) = OP22(:, eidP) + 0.5*K11*Dx(eidP, :)'*diag(-1*nx)*FacialMass2d + 0.5*K13*Dz(eidP, :)'*diag(-1*nx)*FacialMass2d + ...
                    0.5*K22*Dy(eidP,:)'*diag(-1*ny)*FacialMass2d + 0.5*K23*Dz(eidP, :)'*diag(-1*ny)*FacialMass2d;                
    OP21(:, eidP) = OP21(:, eidP) + 0.5*K11*Dx(eidM, :)'*diag(-1*nx)*FacialMass2d + 0.5*K13*Dz(eidM, :)'*diag(-1*nx)*FacialMass2d + ...
                    0.5*K22*Dy(eidM,:)'*diag(-1*ny)*FacialMass2d + 0.5*K23*Dz(eidM, :)'*diag(-1*ny)*FacialMass2d;             
    %> For term $\left {\nabla_y q_h\right }[v]_x$
    OP11(eidM, :) = OP11(eidM, :) + 0.5 * K11 * diag(nx) * FacialMass2d * Dx(eidM, :) + 0.5 * K13*diag(nx)*FacialMass2d*Dz(eidM, :)+...
                    0.5*K22*diag(ny)*FacialMass2d*Dy(eidM,:) + 0.5*K23*diag(ny)*FacialMass2d*Dz(eidM, :);
    OP12(eidP,:) = OP12(eidP, :) + 0.5 * K11 * diag(-1*nx)* FacialMass2d * Dx(eidM, :) + 0.5 * K13 * diag(-1*nx)* FacialMass2d * Dz(eidM, :) + ...
                    0.5*K22*diag(-1*ny)*FacialMass2d*Dy(eidM,:) + 0.5 * K23 * diag(-1*ny)*FacialMass2d*Dz(eidM,:);
    OP22(eidP, :) = OP22(eidP, :) + 0.5 * K11 * diag(-1*nx) * FacialMass2d * Dx(eidP, :) + 0.5 * K13*diag(-1*nx)*FacialMass2d*Dz(eidP, :)+...
                    0.5*K22*diag(-1*ny)*FacialMass2d*Dy(eidP,:) + 0.5*K23*diag(-1*ny)*FacialMass2d*Dz(eidP, :);
    OP21(eidM, :) = OP21(eidM, :) + 0.5 * K11 * diag(nx) * FacialMass2d *Dx(eidP,:) + 0.5*K13*diag(nx)*FacialMass2d*Dz(eidP,:) + ...
                    0.5*K22*diag(ny)*FacialMass2d*Dy(eidP,:) + 0.5*K23*diag(ny)*FacialMass2d*Dz(eidP,:);
    %> For term $-\tau [q_h]_x[v]_x$  check to here
    OP11(eidM, eidM) = OP11(eidM, eidM) - Tau(face) * FacialMass2d;
    OP12(eidP, eidM) = OP12(eidP, eidM) - Tau(face) * (-1) * FacialMass2d;
    OP22(eidP, eidP) = OP22(eidP, eidP) - Tau(face) * FacialMass2d;
    OP21(eidM, eidP) = OP21(eidM, eidP) - Tau(face) * (-1) * FacialMass2d;

    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OP11;
    obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OP12;
    obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = ...
        obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        ElementMassMatrix\OP22;
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        ElementMassMatrix\OP21;
end

fm = repmat( (obj.meshUnion.BottomEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.BottomEdge.FToE(1,:)))', 1, 1 );
Tau = bsxfun(@times,  (fm(:) )',...
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
    %> For term $[q_h]_y\left {\nabla_x v\right }$
    OP11(:, eidM) = OP11(:, eidM) + 0.5*K13*Dx(eidM, :)'*diag(nz)*FacialMass2d + 0.5 * K23*Dy(eidM, :)'*diag(nz)*FacialMass2d + ...
                 0.5 * K33*Dz(eidM,:)'*diag(nz)*FacialMass2d; % Horizontal direction not consider here
                
    OP12(:, eidM) = OP12(:, eidM) + 0.5*K13*Dx(eidP,:)'*diag(nz)*FacialMass2d + 0.5 * K23*Dy(eidP, :)'*diag(nz)*FacialMass2d + ...
                 0.5 * K33*Dz(eidP,:)'*diag(nz)*FacialMass2d;
    OP22(:, eidP) = OP22(:, eidP) + 0.5*K13*Dx(eidP, :)'*diag(-1*nz)*FacialMass2d + 0.5 * K23*Dy(eidP, :)'*diag(-1*nz)*FacialMass2d + ...
                 0.5 * K33*Dz(eidP,:)'*diag(-1*nz)*FacialMass2d;
    OP21(:, eidP) = OP21(:, eidP) + 0.5*K13*Dx(eidM, :)'*diag(-1*nz)*FacialMass2d + 0.5 * K23*Dy(eidM, :)'*diag(-1*nz)*FacialMass2d + ...
                 0.5 * K33*Dz(eidM,:)'*diag(-1*nz)*FacialMass2d; 
    %> For term $\left {\nabla_y q_h\right }[v]_x$
    OP11(eidM, :) = OP11(eidM, :) + 0.5*K13*diag(nz)*FacialMass2d*Dx(eidM, :) + 0.5*K23*diag(nz)*FacialMass2d*Dy(eidM, :) + ...
        0.5 * K33 * diag(nz) * FacialMass2d*Dz(eidM, :);
    OP12(eidP, :) = OP12(eidP, :) + 0.5*K13*diag(-1*nz)*FacialMass2d*Dx(eidM, :) + 0.5*K23*diag(-1*nz)*FacialMass2d*Dy(eidM, :) + ...
        0.5 * K33 * diag(-1*nz) * FacialMass2d*Dz(eidM, :);    
    OP22(eidP, :) = OP22(eidP, :) + 0.5*K13*diag(-1*nz)*FacialMass2d*Dx(eidP, :) + 0.5*K23*diag(-1*nz)*FacialMass2d*Dy(eidP, :) + ...
        0.5 * K33 * diag(-1*nz) * FacialMass2d*Dz(eidP, :);
    OP21(eidM, :) = OP21(eidM, :) + 0.5*K13*diag(nz)*FacialMass2d*Dx(eidP, :) + 0.5*K23*diag(nz)*FacialMass2d*Dy(eidP, :) + ...
        0.5 * K33 * diag(nz) * FacialMass2d*Dz(eidP, :); 
    %> For term $-\tau [q_h]_x[v]_x$  check to here
    OP11(eidM, eidM) = OP11(eidM, eidM) - Tau(face) * FacialMass2d;
    OP12(eidP, eidM) = OP12(eidP, eidM) - Tau(face) * (-1) * FacialMass2d;
    OP22(eidP, eidP) = OP22(eidP, eidP) - Tau(face) * FacialMass2d;
    OP21(eidM, eidP) = OP21(eidM, eidP) - Tau(face) * (-1) * FacialMass2d;
    
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OP11;
    obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OP12;
    obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = ...
        obj.StiffMatrix((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        ElementMassMatrix\OP22;
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        ElementMassMatrix\OP21;
end

fm = repmat( (obj.meshUnion.SurfaceBoundaryEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.SurfaceBoundaryEdge.FToE(1,:)))', 1, 1 );
Tau = bsxfun(@times,  (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    3 )/3 * obj.meshUnion(1).cell.Nface/2);
for face = 1:SurfaceBoundaryEdge.Ne
    OP11 = zeros(Np);
    ele = SurfaceBoundaryEdge.FToE(2*(face-1)+1);
    eidM = SurfaceBoundaryEdge.FToN1(:,face);
    nz = SurfaceBoundaryEdge.nz(:,face);
    Js = SurfaceBoundaryEdge.Js(:,face);
    FacialMass2d = diag(Js)*SurfaceBoundaryEdge.M;
    
    OP11(:, eidM) = OP11(:, eidM) + K13*Dx(eidM, :)'*diag(nz)*FacialMass2d + K23*Dy(eidM, :)'*diag(nz)*FacialMass2d + ...
                    K33*Dz(eidM,:)'*diag(nz)*FacialMass2d; % Horizontal direction not consider here
                
%   %> For term $\nabla_y q_h v n_xd\boldsymbol{x}$
    
    OP11(eidM, :) = OP11(eidM, :) + K13 * diag(nz) * FacialMass2d * Dx(eidM, :) + K23*diag(nz)*FacialMass2d*Dy(eidM, :)+...
                    K33*diag(nz)*FacialMass2d*Dz(eidM,:);
    %> For term $-\tau n_x^2q_h v$
    OP11(eidM,eidM) = OP11(eidM,eidM) - Tau(face)*FacialMass2d;
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OP11;
    
end
obj.StiffMatrix = sparse(obj.StiffMatrix);
end