function AssembleGlobalStiffMatrix(obj)
%%%%%%%Not implemented %%%%%%%%%%%%%%%%%%%%%%%%
% [ TempRHS, TempStiffMatrix ] = mxAssembleGlobalStiffMatrixWithBCsImposed(obj.NonhydrostaticSolver.MSPNPX + obj.NonhydrostaticSolver.MSPNPY, ...
%    obj.RHS, obj.BoundaryDirichletData, obj.NonhydrostaticSolver.BoundaryEdge, obj.NonhydrostaticSolver.cell, obj.NonhydrostaticSolver.mesh, ...
%    obj.NonhydrostaticSolver.SurfaceBoundaryEdge, obj.NonhydrostaticSolver.BottomBoundaryEdge, obj.SurfaceDirichletData, ...
%    obj.BottomBoundaryDirichletData);

% For this case we only test where the matlab version and C version is
% equal when periodic boundary condition imposed at the lateral boundary, Newmann boundary
% imposed at the bottom boundary and Dirichlet boundary condition imposed on the
% surface boundary
TempStiffMatrix = obj.NonhydrostaticSolver.MSPNPX + obj.NonhydrostaticSolver.MSPNPY;

AssembleMixedOrderTermGlobalStiffMatrix(obj);

disp("=======================The maximum difference is====================");
disp(max(max(TempStiffMatrix - obj.StiffMatrix)));
disp("=======================The minimum difference is====================");
disp(min(min(TempStiffMatrix - obj.StiffMatrix)));
end

function AssembleMixedOrderTermGlobalStiffMatrix(obj)
%> For term $\frac{\partial}{\partial x}\left (\frac{\partial u}{\partial y}\right )$
mesh = obj.meshUnion;
InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;
BottomEdge = mesh.BottomEdge;
SurfaceBoundaryEdge = mesh.SurfaceBoundaryEdge;
BottomBoundaryEdge = mesh.BottomBoundaryEdge;
cell = obj.meshUnion.cell;
K = mesh.K;  Np = cell.Np;
Nfp = InnerEdge.Nfp;
obj.StiffMatrix = zeros(K*Np);
%For this test, the mesh is uniform
Dx = diag(mesh.rx(:,1))*mesh.cell.Dr +  diag(mesh.sx(:,1))*mesh.cell.Ds +  diag(mesh.tx(:,1))*mesh.cell.Dt;
Dy = diag(mesh.ry(:,1))*mesh.cell.Dr +  diag(mesh.sy(:,1))*mesh.cell.Ds +  diag(mesh.ty(:,1))*mesh.cell.Dt;
Dz = diag(mesh.rz(:,1))*mesh.cell.Dr +  diag(mesh.sz(:,1))*mesh.cell.Ds +  diag(mesh.tz(:,1))*mesh.cell.Dt;
ElementMassMatrix = diag(mesh.J(:,1))*cell.M;
%> For term $-\int_{\Omega_e}\nabla_x v\nabla_y q_hd\boldsymbol {x}$
OP11 = -1 * ( ElementMassMatrix\(Dx'*ElementMassMatrix*Dz) ) - ...
    ( ElementMassMatrix\(Dy'*ElementMassMatrix*Dz) ) ;  %checked
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
    nz = InnerEdge.nz(:,face);
    Js = InnerEdge.Js(:,face);
    FacialMass2d = diag(Js)*InnerEdge.M;
    %> For term $[q_h]_y\left {\nabla_x v\right }$
    OP11(:, eidM) = OP11(:, eidM) + 0.5*Dx(eidM,:)'*diag(nz)*FacialMass2d + 0.5*Dy(eidM,:)'*diag(nz)*FacialMass2d;   %checked
    OP12(:, eidM) = OP12(:, eidM) + 0.5*Dx(eidP,:)'*diag(nz)*FacialMass2d + 0.5*Dy(eidP,:)'*diag(nz)*FacialMass2d;   %checked
    OP22(:, eidP) = OP22(:, eidP) + 0.5*Dx(eidP,:)'*diag(-1*nz)*FacialMass2d + 0.5*Dy(eidP,:)'*diag(-1*nz)*FacialMass2d;
    OP21(:, eidP) = OP21(:, eidP) + 0.5*Dx(eidM,:)'*diag(-1*nz)*FacialMass2d + 0.5*Dy(eidM,:)'*diag(-1*nz)*FacialMass2d;
    %> For term $\left {\nabla_y q_h\right }[v]_x$
    OP11(eidM, :) = OP11(eidM, :) + 0.5 * diag(nx)*FacialMass2d*Dz(eidM,:) + 0.5 * diag(ny)*FacialMass2d*Dz(eidM,:);
    OP12(eidP, :) = OP12(eidP, :) + 0.5 * diag(-1*nx)*FacialMass2d*Dz(eidM,:) + 0.5 * diag(-1*ny)*FacialMass2d*Dz(eidM,:);
    OP22(eidP, :) = OP22(eidP, :) + 0.5 * diag(-1*nx)*FacialMass2d*Dz(eidP,:) + 0.5 * diag(-1*ny)*FacialMass2d*Dz(eidP,:);
    OP21(eidM, :) = OP21(eidM, :) + 0.5 * diag(nx)*FacialMass2d*Dz(eidP,:) + 0.5 * diag(ny)*FacialMass2d*Dz(eidP,:);
    %> For term $-\tau [q_h]_x[v]_x$  check to here
    %     OP11(eidM, eidM) = OP11(eidM, eidM) - Tau(face)* diag(nx.^2) * FacialMass2d;
    %     OP12(eidP, eidM) = OP12(eidP, eidM) - Tau(face)* diag(-1*nx.^2) * FacialMass2d;
    %     OP22(eidP, eidP) = OP22(eidP, eidP) - Tau(face)* diag(nx.^2) * FacialMass2d;
    %     OP21(eidM, eidP) = OP21(eidM, eidP) - Tau(face)* diag(-1*nx.^2) * FacialMass2d;
    OP11(eidM, eidM) = OP11(eidM, eidM) - 2 * Tau(face)* 1 * FacialMass2d;
    OP12(eidP, eidM) = OP12(eidP, eidM) - 2 * Tau(face)* -1 * FacialMass2d;
    OP22(eidP, eidP) = OP22(eidP, eidP) - 2 * Tau(face)* 1 * FacialMass2d;
    OP21(eidM, eidP) = OP21(eidM, eidP) - 2 * Tau(face)* -1 * FacialMass2d;
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
    nx = BottomEdge.nx(:,face);
    ny = BottomEdge.ny(:,face);
    nz = BottomEdge.nz(:,face);
    Js = BottomEdge.Js(:,face);
    FacialMass2d = diag(Js)*BottomEdge.M;
    %> For term $[q_h]_y\left {\nabla_x v\right }$
    OP11(:, eidM) = OP11(:, eidM) + 0.5*Dx(eidM,:)'*diag(nz)*FacialMass2d + 0.5*Dy(eidM,:)'*diag(nz)*FacialMass2d;   %checked
    OP12(:, eidM) = OP12(:, eidM) + 0.5*Dx(eidP,:)'*diag(nz)*FacialMass2d + 0.5*Dy(eidP,:)'*diag(nz)*FacialMass2d;   %checked
    OP22(:, eidP) = OP22(:, eidP) + 0.5*Dx(eidP,:)'*diag(-1*nz)*FacialMass2d + 0.5*Dy(eidP,:)'*diag(-1*nz)*FacialMass2d;
    OP21(:, eidP) = OP21(:, eidP) + 0.5*Dx(eidM,:)'*diag(-1*nz)*FacialMass2d + 0.5*Dy(eidM,:)'*diag(-1*nz)*FacialMass2d;
    %> For term $\left {\nabla_y q_h\right }[v]_x$
    OP11(eidM, :) = OP11(eidM, :) + 0.5 * diag(nx)*FacialMass2d*Dz(eidM,:) + 0.5 * diag(ny)*FacialMass2d*Dz(eidM,:);
    OP12(eidP, :) = OP12(eidP, :) + 0.5 * diag(-1*nx)*FacialMass2d*Dz(eidM,:) + 0.5 * diag(-1*ny)*FacialMass2d*Dz(eidM,:);
    OP22(eidP, :) = OP22(eidP, :) + 0.5 * diag(-1*nx)*FacialMass2d*Dz(eidP,:) + 0.5 * diag(-1*ny)*FacialMass2d*Dz(eidP,:);
    OP21(eidM, :) = OP21(eidM, :) + 0.5 * diag(nx)*FacialMass2d*Dz(eidP,:) + 0.5 * diag(ny)*FacialMass2d*Dz(eidP,:);
    %> For term $-\tau [q_h]_x[v]_x$  check to here
    %     OP11(eidM, eidM) = OP11(eidM, eidM) - Tau(face)* diag(nx.^2) * FacialMass2d;
    %     OP12(eidP, eidM) = OP12(eidP, eidM) - Tau(face)* diag(-1*nx.^2) * FacialMass2d;
    %     OP22(eidP, eidP) = OP22(eidP, eidP) - Tau(face)* diag(nx.^2) * FacialMass2d;
    %     OP21(eidM, eidP) = OP21(eidM, eidP) - Tau(face)* diag(-1*nx.^2) * FacialMass2d;
    OP11(eidM, eidM) = OP11(eidM, eidM) - 2 * Tau(face)* 1 * FacialMass2d;
    OP12(eidP, eidM) = OP12(eidP, eidM) - 2 * Tau(face)* -1 * FacialMass2d;
    OP22(eidP, eidP) = OP22(eidP, eidP) - 2 * Tau(face)* 1 * FacialMass2d;
    OP21(eidM, eidP) = OP21(eidM, eidP) - 2 * Tau(face)* -1 * FacialMass2d;
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
    nx = SurfaceBoundaryEdge.nx(:,face);
    ny = SurfaceBoundaryEdge.ny(:,face);
    nz = SurfaceBoundaryEdge.nz(:,face);
    Js = SurfaceBoundaryEdge.Js(:,face);
    FacialMass2d = diag(Js)*SurfaceBoundaryEdge.M;
    % At present, at surface, we only impose Dirichlet boundary condition
    %> For term $q_D\nabla_xv n_y d\boldsymbol{x}$
    TempData = sum( Dx(eidM,:)' * FacialMass2d * diag(nz) * diag(obj.SurfaceDirichletData((face-1)*Nfp+1:face*Nfp)) , 2) + ...
        sum( Dy(eidM,:)' * FacialMass2d * diag(nz) * diag(obj.SurfaceDirichletData((face-1)*Nfp+1:face*Nfp)) , 2);
    obj.RHS(:, ele) = obj.RHS(:, ele) + ( ElementMassMatrix \ TempData );
    %> For term $q_h\nabla_x vn_yd\boldsymbol{x}$
    OP11(:,eidM) = OP11(:,eidM) + Dx(eidM,:)'*diag(nz)*FacialMass2d + Dy(eidM,:)'*diag(nz)*FacialMass2d;
    %> For term $\nabla_y q_h v n_xd\boldsymbol{x}$
    OP11(eidM,:) = OP11(eidM,:) + diag(nx)*FacialMass2d * Dz(eidM,:) + diag(ny)*FacialMass2d * Dz(eidM,:);
    %> For term $-\tau n_x^2q_h v$
    %             OP11(eidM,eidM) = OP11(eidM,eidM) - diag(nx.^2)*Tau(face)*FacialMass2d;
    OP11(eidM,eidM) = OP11(eidM,eidM) - 2 * Tau(face)*FacialMass2d;
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OP11;
    %> For term $-\tau n_x^2 q_D v d\boldsymbol{x}$
    TempData = zeros(Np,1);
    %             TempData(eidM) = sum( Tau(face) * FacialMass2d * diag(nx.^2)*...
    %                 diag(obj.DirichletData((face-1)*Nfp+1:face*Nfp)), 2 );
    TempData(eidM) = sum( 2 * Tau(face) * FacialMass2d*...
        diag(obj.SurfaceDirichletData((face-1)*Nfp+1:face*Nfp)), 2 );
    obj.RHS(:,ele) = obj.RHS(:, ele) - ( ElementMassMatrix \ TempData );
    
end

fm = repmat( (obj.meshUnion.BoundaryEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.BoundaryEdge.FToE(1,:)))', 1, 1 );
Tau = bsxfun(@times,  (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    3 )/3 * obj.meshUnion(1).cell.Nface/2);

for face = 1:BoundaryEdge.Ne
    OP11 = zeros(Np);
    ele = BoundaryEdge.FToE(2*(face-1)+1);
    eidM = BoundaryEdge.FToN1(:,face);
    nx = BoundaryEdge.nx(:,face);
    ny = BoundaryEdge.ny(:,face);
    nz = BoundaryEdge.nz(:,face);
    Js = BoundaryEdge.Js(:,face);
    FacialMass2d = diag(Js)*BoundaryEdge.M;
    % At present, at surface, we only impose Dirichlet boundary condition
    %> For term $q_D\nabla_xv n_y d\boldsymbol{x}$
    TempData = sum( Dx(eidM,:)' * FacialMass2d * diag(nz) * diag(obj.BoundaryDirichletData((face-1)*Nfp+1:face*Nfp)) , 2) + ...
        sum( Dy(eidM,:)' * FacialMass2d * diag(nz) * diag(obj.BoundaryDirichletData((face-1)*Nfp+1:face*Nfp)) , 2);
    obj.RHS(:, ele) = obj.RHS(:, ele) + ( ElementMassMatrix \ TempData );
    %> For term $q_h\nabla_x vn_yd\boldsymbol{x}$
    OP11(:,eidM) = OP11(:,eidM) + Dx(eidM,:)'*diag(nz)*FacialMass2d + Dy(eidM,:)'*diag(nz)*FacialMass2d;
    %> For term $\nabla_y q_h v n_xd\boldsymbol{x}$
    OP11(eidM,:) = OP11(eidM,:) + diag(nx)*FacialMass2d * Dz(eidM,:) + diag(ny)*FacialMass2d * Dz(eidM,:);
    %> For term $-\tau n_x^2q_h v$
    %             OP11(eidM,eidM) = OP11(eidM,eidM) - diag(nx.^2)*Tau(face)*FacialMass2d;
    OP11(eidM,eidM) = OP11(eidM,eidM) - 2 * Tau(face)*FacialMass2d;
    obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OP11;
    %> For term $-\tau n_x^2 q_D v d\boldsymbol{x}$
    TempData = zeros(Np,1);
    %             TempData(eidM) = sum( Tau(face) * FacialMass2d * diag(nx.^2)*...
    %                 diag(obj.DirichletData((face-1)*Nfp+1:face*Nfp)), 2 );
    TempData(eidM) = sum( 2 * Tau(face) * FacialMass2d*...
        diag(obj.BoundaryDirichletData((face-1)*Nfp+1:face*Nfp)), 2 );
    obj.RHS(:,ele) = obj.RHS(:, ele) - ( ElementMassMatrix \ TempData );
    
end

fm = repmat( (obj.meshUnion.BottomBoundaryEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.BottomBoundaryEdge.FToE(1,:)))', 1, 1 );
Tau = bsxfun(@times,  (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    3 )/3 * obj.meshUnion(1).cell.Nface/2);

% for face = 1:BottomBoundaryEdge.Ne
%     OP11 = zeros(Np);
%     ele = BottomBoundaryEdge.FToE(2*(face-1)+1);
%     eidM = BottomBoundaryEdge.FToN1(:,face);
%     nx = BottomBoundaryEdge.nx(:,face);
%     ny = BottomBoundaryEdge.ny(:,face);
%     nz = BottomBoundaryEdge.nz(:,face);
%     Js = BottomBoundaryEdge.Js(:,face);
%     FacialMass2d = diag(Js)*BottomBoundaryEdge.M;
%     % At present, at surface, we only impose Dirichlet boundary condition
%     %> For term $q_D\nabla_xv n_y d\boldsymbol{x}$
%     TempData = sum( Dx(eidM,:)' * FacialMass2d * diag(nz) * diag(obj.BottomBoundaryDirichletData((face-1)*Nfp+1:face*Nfp)) , 2) + ...
%         sum( Dy(eidM,:)' * FacialMass2d * diag(nz) * diag(obj.BottomBoundaryDirichletData((face-1)*Nfp+1:face*Nfp)) , 2);
%     obj.RHS(:, ele) = obj.RHS(:, ele) + ( ElementMassMatrix \ TempData );
%     %> For term $q_h\nabla_x vn_yd\boldsymbol{x}$
%     OP11(:,eidM) = OP11(:,eidM) + Dx(eidM,:)'*diag(nz)*FacialMass2d + Dy(eidM,:)'*diag(nz)*FacialMass2d;
%     %> For term $\nabla_y q_h v n_xd\boldsymbol{x}$
%     OP11(eidM,:) = OP11(eidM,:) + diag(nx)*FacialMass2d * Dz(eidM,:) + diag(ny)*FacialMass2d * Dz(eidM,:);
%     %> For term $-\tau n_x^2q_h v$
%     %             OP11(eidM,eidM) = OP11(eidM,eidM) - diag(nx.^2)*Tau(face)*FacialMass2d;
%     OP11(eidM,eidM) = OP11(eidM,eidM) - 2 * Tau(face)*FacialMass2d;
%     obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
%         ElementMassMatrix\OP11;
%     %> For term $-\tau n_x^2 q_D v d\boldsymbol{x}$
%     TempData = zeros(Np,1);
%     %             TempData(eidM) = sum( Tau(face) * FacialMass2d * diag(nx.^2)*...
%     %                 diag(obj.DirichletData((face-1)*Nfp+1:face*Nfp)), 2 );
%     TempData(eidM) = sum( 2 * Tau(face) * FacialMass2d*...
%         diag(obj.BottomBoundaryDirichletData((face-1)*Nfp+1:face*Nfp)), 2 );
%     obj.RHS(:,ele) = obj.RHS(:, ele) - ( ElementMassMatrix \ TempData );
%     
% end

obj.StiffMatrix = sparse(obj.StiffMatrix);
end










