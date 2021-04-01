function AssembleGlobalStiffMatrix(obj)

AssembleMixedOrderTermGlobalStiffMatrix(obj);

end

function AssembleMixedOrderTermGlobalStiffMatrix(obj)
%> For term $\frac{\partial}{\partial x}\left (\frac{\partial u}{\partial y}\right ) + \frac{\partial}{\partial y}\left (\frac{\partial u}{\partial x}\right )$
mesh = obj.meshUnion;
InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;
cell = obj.meshUnion.cell;
K = mesh.K;  Np = cell.Np;
Nfp = InnerEdge.Nfp;
obj.StiffMatrix = zeros(K*Np);
%For this test, the mesh is uniform
Dx = diag(mesh.rx(:,1))*mesh.cell.Dr +  diag(mesh.sx(:,1))*mesh.cell.Ds;
Dy = diag(mesh.ry(:,1))*mesh.cell.Dr +  diag(mesh.sy(:,1))*mesh.cell.Ds;
ElementMassMatrix = diag(mesh.J(:,1))*cell.M;
%> For term $-\int_{\Omega_e}\nabla_x v\nabla_y q_hd\boldsymbol {x}$
OP11 = -1 * (ElementMassMatrix\(Dx'*ElementMassMatrix*Dy) ) -1 * (ElementMassMatrix\(Dy'*ElementMassMatrix*Dx))...
    -1 * (ElementMassMatrix\(Dx' * ElementMassMatrix * Dx + Dy' * ElementMassMatrix * Dy));  %checked
for i = 1:K
    obj.StiffMatrix((i-1)*Np+1:i*Np,(i-1)*Np+1:i*Np) = OP11;
end
fm = repmat( (obj.meshUnion.InnerEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.InnerEdge.FToE(1,:)))', 1, 1 );
Tau12 = bsxfun(@times,  obj.D12 * (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    2 )/2 * obj.meshUnion(1).cell.Nface/2);
Tau21 = bsxfun(@times,  obj.D21*(fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    2 )/2 * obj.meshUnion(1).cell.Nface/2);
Tau1122 = bsxfun(@times,  obj.D11*(fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    2 )/2 * obj.meshUnion(1).cell.Nface/2);

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
    OP11(:, eidM) = OP11(:, eidM) + 0.5*Dx(eidM,:)'*diag(ny)*FacialMass2d+ 0.5*Dy(eidM,:)'*diag(nx)*FacialMass2d + 0.5 * ( Dx(eidM,:) )' * ( FacialMass2d ) * diag(nx) + 0.5 * ( Dy(eidM,:) )' * ( FacialMass2d ) * diag(ny);   %checked
    OP12(:, eidM) = OP12(:, eidM) + 0.5*Dx(eidP,:)'*diag(ny)*FacialMass2d+ 0.5*Dy(eidP,:)'*diag(nx)*FacialMass2d + 0.5 * ( Dx(eidP,:) )' * ( FacialMass2d ) * diag(nx) + 0.5 * ( Dy(eidP,:) )' * ( FacialMass2d ) * diag(ny);   %checked
    OP22(:, eidP) = OP22(:, eidP) + 0.5*Dx(eidP,:)'*diag(-1*ny)*FacialMass2d+ 0.5*Dy(eidP,:)'*diag(-1*nx)*FacialMass2d + 0.5 * ( Dx(eidP,:) )' * ( FacialMass2d ) * diag(-1*nx) + 0.5 * ( Dy(eidP,:) )' * ( FacialMass2d ) * diag(-1*ny);
    OP21(:, eidP) = OP21(:, eidP) + 0.5*Dx(eidM,:)'*diag(-1*ny)*FacialMass2d+ 0.5*Dy(eidM,:)'*diag(-1*nx)*FacialMass2d + 0.5 * ( Dx(eidM,:) )' * ( FacialMass2d ) * diag(-1*nx) + 0.5 * ( Dy(eidM,:) )' * ( FacialMass2d ) * diag(-1*ny);
    %> For term $\left {\nabla_y q_h\right }[v]_x$
    OP11(eidM, :) = OP11(eidM, :) + 0.5 * diag(nx)*FacialMass2d*Dy(eidM,:)+ 0.5 * diag(ny)*FacialMass2d*Dx(eidM,:) + 0.5 * FacialMass2d * diag(nx) * Dx(eidM,:) + 0.5 * FacialMass2d * diag(ny) * Dy(eidM,:);
    OP12(eidP, :) = OP12(eidP, :) + 0.5 * diag(-1*nx)*FacialMass2d*Dy(eidM,:)+ 0.5 * diag(-1*ny)*FacialMass2d*Dx(eidM,:) + 0.5 * FacialMass2d * diag(-1*nx) * Dx(eidM,:) + 0.5 * FacialMass2d * diag(-1*ny) * Dy(eidM,:);
    OP22(eidP, :) = OP22(eidP, :) + 0.5 * diag(-1*nx)*FacialMass2d*Dy(eidP,:)+ 0.5 * diag(-1*ny)*FacialMass2d*Dx(eidP,:) + 0.5 * FacialMass2d * diag(-1*nx) * Dx(eidP,:) + 0.5 * FacialMass2d * diag(-1*ny) * Dy(eidP,:);
    OP21(eidM, :) = OP21(eidM, :) + 0.5 * diag(nx)*FacialMass2d*Dy(eidP,:)+ 0.5 * diag(ny)*FacialMass2d*Dx(eidP,:) + 0.5 * FacialMass2d * diag(nx) * Dx(eidP,:) + 0.5 * FacialMass2d * diag(ny) * Dy(eidP,:);
    %> For term $-\tau [q_h]_x[v]_x$  check to here
    OP11(eidM, eidM) = OP11(eidM, eidM) - Tau12(face)* 1 * FacialMass2d;
    OP12(eidP, eidM) = OP12(eidP, eidM) - Tau12(face)* -1 * FacialMass2d;
    OP22(eidP, eidP) = OP22(eidP, eidP) - Tau12(face)* 1 * FacialMass2d;
    OP21(eidM, eidP) = OP21(eidM, eidP) - Tau12(face)* -1 * FacialMass2d;
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

fm = repmat( (obj.meshUnion.BoundaryEdge.LAV./obj.meshUnion.LAV(obj.meshUnion.BoundaryEdge.FToE(1,:)))', 1, 1 );
Tau12 = bsxfun(@times,  obj.D12 * (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    2 )/2 * obj.meshUnion(1).cell.Nface/2);
Tau21 = bsxfun(@times,  obj.D21 * (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    2 )/2 * obj.meshUnion(1).cell.Nface/2);
Tau1122 = bsxfun(@times,  obj.D11 * (fm(:) )',...
    ( obj.meshUnion(1).cell.N + 1 )*(obj.meshUnion(1).cell.N + ...
    2 )/2 * obj.meshUnion(1).cell.Nface/2);
for face = 1:BoundaryEdge.Ne
    OP11 = zeros(Np);
    ele = BoundaryEdge.FToE(2*(face-1)+1);
    eidM = BoundaryEdge.FToN1(:,face);
    nx = BoundaryEdge.nx(:,face);
    ny = BoundaryEdge.ny(:,face);
    Js = BoundaryEdge.Js(:,face);
    FacialMass2d = diag(Js)*BoundaryEdge.M;
    switch(BoundaryEdge.ftype(face))
        case enumBoundaryCondition.Newmann
            TempData = zeros(Np,1);
            TempData(eidM) = FacialMass2d * obj.NewmannData(:,face);
            obj.RHS((ele-1)*Np+1:ele*Np) = obj.RHS((ele-1)*Np+1:ele*Np) - ( ElementMassMatrix \ TempData )';
        case enumBoundaryCondition.Dirichlet
            %> For term $q_D\nabla_xv n_y d\boldsymbol{x}$
            TempData = sum( obj.D12 * Dx(eidM,:)' * FacialMass2d * diag(ny) * diag(obj.DirichletData((face-1)*Nfp+1:face*Nfp))  + ...
                obj.D21 * Dy(eidM,:)' * FacialMass2d * diag(nx) * diag(obj.DirichletData((face-1)*Nfp+1:face*Nfp)), 2);
            obj.RHS((ele-1)*Np+1:ele*Np) = obj.RHS((ele-1)*Np+1:ele*Np) + ( ElementMassMatrix \ TempData )';
            %> For term $q_h\nabla_x vn_yd\boldsymbol{x}$
            OP11(:,eidM) = OP11(:,eidM) + Dx(eidM,:)'*diag(ny)*FacialMass2d  + Dy(eidM,:)'*diag(nx)*FacialMass2d +  ( Dx(eidM,:) )'* diag(nx) * ( FacialMass2d ) +  ( Dy(eidM,:) )'* diag(ny) * ( FacialMass2d );
            %> For term $\nabla_y q_h v n_xd\boldsymbol{x}$
            OP11(eidM,:) = OP11(eidM,:) + diag(nx)*FacialMass2d * Dy(eidM,:)  + diag(ny)*FacialMass2d * Dx(eidM,:) +  FacialMass2d * diag(nx) * Dx(eidM,:) +  FacialMass2d * diag(ny) * Dy(eidM,:);
            %> For term $-\tau n_x^2q_h v$
            OP11(eidM,eidM) = OP11(eidM,eidM) - diag(nx.^2)*Tau12(face)*FacialMass2d - diag(ny.^2)*Tau21(face)*FacialMass2d -  Tau1122(face) * FacialMass2d;
            obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = obj.StiffMatrix((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
                ElementMassMatrix\OP11;
            %> For term $-\tau n_x^2 q_D v d\boldsymbol{x}$
            TempData = zeros(Np,1);
            TempData(eidM) = sum( Tau12(face) * FacialMass2d * diag(nx.^2)*...
                diag(obj.DirichletData((face-1)*Nfp+1:face*Nfp)) + Tau21(face) * FacialMass2d * diag(ny.^2)*...
                diag(obj.DirichletData((face-1)*Nfp+1:face*Nfp)), 2 );
            obj.RHS((ele-1)*Np+1:ele*Np) = obj.RHS((ele-1)*Np+1:ele*Np) - ( ElementMassMatrix \ TempData )';
    end
end
end