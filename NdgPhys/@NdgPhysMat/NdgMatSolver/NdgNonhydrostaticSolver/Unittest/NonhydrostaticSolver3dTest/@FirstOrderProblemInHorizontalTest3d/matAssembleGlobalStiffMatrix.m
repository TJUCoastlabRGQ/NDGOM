function [ StiffMatrixInX, StiffMatrixInY, LBRHS ] = matAssembleGlobalStiffMatrix(obj)
LBRHS = zeros(1);
mesh = obj.meshUnion;
InnerEdge = mesh.InnerEdge;
cell = mesh.cell;
K = mesh.K;
Np = cell.Np;

StiffMatrixInX     = zeros(K*Np);
StiffMatrixInY     = zeros(K*Np);

%For this test, the mesh is uniform
Dx = diag(mesh.rx(:,1))*mesh.cell.Dr +  diag(mesh.sx(:,1))*mesh.cell.Ds +  diag(mesh.tx(:,1))*mesh.cell.Dt;
Dy = diag(mesh.ry(:,1))*mesh.cell.Dr +  diag(mesh.sy(:,1))*mesh.cell.Ds +  diag(mesh.ty(:,1))*mesh.cell.Dt;
ElementMassMatrix = diag(mesh.J(:,1))*cell.M;

for i = 1:K
    StiffMatrixInX((i-1)*Np+1:i*Np,(i-1)*Np+1:i*Np) = Dx;
    StiffMatrixInY((i-1)*Np+1:i*Np,(i-1)*Np+1:i*Np) = Dy;
end

for face = 1:InnerEdge.Ne
    OPX11 = zeros(Np);
    OPX12 = zeros(Np);
    OPX22 = zeros(Np);
    OPX21 = zeros(Np);
    
    OPY11 = zeros(Np);
    OPY12 = zeros(Np);
    OPY22 = zeros(Np);
    OPY21 = zeros(Np);
    
    ele = InnerEdge.FToE(2*(face-1)+1);
    adjacentEle = InnerEdge.FToE(2*(face-1)+2);
    eidM = InnerEdge.FToN1(:,face);
    eidP = InnerEdge.FToN2(:,face);
    nx = InnerEdge.nx(:,face);
    ny = InnerEdge.ny(:,face);
    Js = InnerEdge.Js(:,face);
    FacialMass2d = diag(Js)*InnerEdge.M;
    
    OPX11(eidM, eidM) = OPX11(eidM, eidM) - 0.5*diag(nx)*FacialMass2d;
    OPY11(eidM, eidM) = OPY11(eidM, eidM) - 0.5*diag(ny)*FacialMass2d;
    
    OPX22(eidP, eidP) = OPX22(eidP, eidP) - 0.5*diag(-1*nx)*FacialMass2d;
    OPY22(eidP, eidP) = OPY22(eidP, eidP) - 0.5*diag(-1*ny)*FacialMass2d;
    
    OPX12(eidP, eidM) = OPX12(eidP, eidM) - 0.5*diag(nx)*FacialMass2d;
    OPY12(eidP, eidM) = OPY12(eidP, eidM) - 0.5*diag(ny)*FacialMass2d;
    
    OPX21(eidM, eidP) = OPX21(eidM, eidP) - 0.5*diag(-1*nx)*FacialMass2d;
    OPY21(eidM, eidP) = OPY21(eidM, eidP) - 0.5*diag(-1*ny)*FacialMass2d;
    
    StiffMatrixInX((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = StiffMatrixInX((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OPX11;
    StiffMatrixInX((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) = StiffMatrixInX((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OPX12;
    StiffMatrixInX((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = ...
        StiffMatrixInX((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        ElementMassMatrix\OPX22;
    StiffMatrixInX((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = StiffMatrixInX((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        ElementMassMatrix\OPX21;
    
    StiffMatrixInY((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) = StiffMatrixInY((ele-1)*Np+1:ele*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OPY11;
    StiffMatrixInY((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) = StiffMatrixInY((adjacentEle-1)*Np+1:adjacentEle*Np,(ele-1)*Np+1:ele*Np) + ...
        ElementMassMatrix\OPY12;
    StiffMatrixInY((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = ...
        StiffMatrixInY((adjacentEle-1)*Np+1:adjacentEle*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        ElementMassMatrix\OPY22;
    StiffMatrixInY((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) = StiffMatrixInY((ele-1)*Np+1:ele*Np,(adjacentEle-1)*Np+1:adjacentEle*Np) + ...
        ElementMassMatrix\OPY21;
end


end