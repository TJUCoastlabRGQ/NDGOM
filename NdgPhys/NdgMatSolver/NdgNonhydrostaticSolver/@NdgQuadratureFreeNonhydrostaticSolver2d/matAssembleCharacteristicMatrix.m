function  [qx, qy, SecondOrderTerm ] = matAssembleCharacteristicMatrix(obj, mesh, ele, edgeType)
%>@brief function used to calculate the term used when assemble the global stiff matrix
%>@detail In this version, the zero boundary condition for nonhydrostatic pressure is imposed at
%> the general clamped boundary, and the zero grad boundary condition at the wall boundary.
%> Also the zero gradient boundary condition for the auxialary variable is imposed at the
%> general clamped boundary, and the zero boundary condition at the wall
%> boundary. This script is used to get the discretization corresponding to
%> $-(\nabla_h s, \nabla_h p_h)_{\Omega}...    
%> +<[p_h],{\nabla_h s}>_{\epsilon_i}...          
%> +<{\nabla_h p_h}-\tau^k[p_h],[s]>_{\epsilon_i}...   
%> -<p_D-p_h,\nabla_h s\cdot \bold{n}>_{\partial \Omega^D}...  
%> +<\bold{n}\cdot \nabla_h p_h,s>_{\partial \Omega^D}...
%> -<\tau^k p_h, s>_{\partial \Omega^D}...    
%> +<\tau^k p_D, s>_{\partial \Omega^D}...
%> + <s, p_N>_{\partial \Omega^N}...
%> =(f,s)_{\Omega}$
%> @param[in] mesh The mesh object
%> @param[in] ele The number of the studied element
%> @param[in] edgeType The edgetype of the studied element
%> @param[out] qx first derivative of nonhydrostatic pressure with respect to x
%> @param[out] qy first derivative of nonhydrostatic pressure with respect to y
%> @param[out] SencondOrderTerm Laplace operator regards to non-hydrostatic pressure
K = mesh.K;  Np = mesh.cell.Np;
SecondOrderTerm = zeros(K*Np, Np);
qx = zeros(K*Np, Np);
qy = zeros(K*Np, Np);

localrows = ((ele-1)*Np+1:ele*Np)';
% cols = ones(1,Np)'*(1:Np);
cols = 1:Np;

Dx = diag(mesh.rx(:,ele))*mesh.cell.Dr +  diag(mesh.sx(:,ele))*mesh.cell.Ds;
Dy = diag(mesh.ry(:,ele))*mesh.cell.Dr +  diag(mesh.sy(:,ele))*mesh.cell.Ds;
ElementMassMatrix = diag(mesh.J(:,ele))*mesh.cell.M;
%> $(\nabla_h s, \nabla_h p_h)_{\Omega}$, this corresponds to the first
%> term in the governing equation, we omit the minus here, and add it at the final part 
OP11 = Dx' * ElementMassMatrix * Dx + Dy' * ElementMassMatrix * Dy;

QX11 = Dx;
QY11 = Dy;

for face = 1:numel(edgeType)
    adjacentEle = mesh.EToE(face, ele);
    adjacentFace = mesh.EToF(face, ele);
%     adjacentrows = ((adjacentEle-1)*Np+1:adjacentEle*Np)'*ones(1,Np);
    adjacentrows = ((adjacentEle-1)*Np+1:adjacentEle*Np)';
    eidM = mesh.cell.Fmask(:,face);
    eidP = mesh.cell.Fmask(:,adjacentFace);
    eidP(1:numel(eidP)) = eidP(numel(eidP):-1:1);
    [nx, ny, Js] = obj.matGetElementFaceNormalDirectionVector( ele, adjacentEle, face,  mesh.BoundaryEdge.FToF, mesh.BoundaryEdge.nx, mesh.BoundaryEdge.ny,...
        mesh.BoundaryEdge.FToE, mesh.InnerEdge.nx, mesh.InnerEdge.ny, mesh.InnerEdge.FToE, mesh.BoundaryEdge.Js, mesh.InnerEdge.Js );
    Dx2 = diag(mesh.rx(:,adjacentEle))*mesh.cell.Dr +  diag(mesh.sx(:,adjacentEle))*mesh.cell.Ds;
    Dy2 = diag(mesh.ry(:,adjacentEle))*mesh.cell.Dr +  diag(mesh.sy(:,adjacentEle))*mesh.cell.Ds;
    
    Dn1 = nx*Dx  + ny*Dy ;
    Dn2 = nx*Dx2 + ny*Dy2;
    
    EdgeMassMatrix = zeros(Np);
    switch(edgeType(face))
        case 2    %For Dirichlet edge boundary
            %> $(\tau^k p_h,s)_{\partial \Omega^D}$, this corresponds to the sixth term in the governning equation, and we omit the minus here
            OP11(eidM, eidM) = OP11(eidM, eidM) +  obj.Tau * Js * mesh.BoundaryEdge.M;
            %> $(s, n\cdot \Nabla p_h)_{\partial \Omega^D}$, this corresponds to the fifth term in the governning equation, and we omit the minus here
            OP11(eidM, :) = OP11(eidM, :) -  Js * mesh.BoundaryEdge.M * Dn1(eidM,:);
            %> $(p_h, \nabla_h s\cdot n)$, this corresponds to the last part in the fourth term in the governning equation, and we omit the minus here
            OP11(:, eidM) = OP11(:, eidM) -  ( Dn1(eidM,:) )' * ( Js * mesh.BoundaryEdge.M );
            
            EdgeMassMatrix(eidM, eidM) = Js * mesh.BoundaryEdge.M;
            QX11 = QX11 - ElementMassMatrix\( nx * EdgeMassMatrix );
            QY11 = QY11 - ElementMassMatrix\( ny * EdgeMassMatrix );
        case 1   %For Newmann edge boundary
            
        otherwise
            %> $(\tau^k [p_h],[s])_{\epsilon_i}$, this corresponds to the last part in the third term in the governning equation, and we omit the minus here
            OP11(eidM, eidM) = OP11(eidM, eidM) +  obj.Tau* Js * mesh.InnerEdge.M;
            %> $<{\nabla_h p_h},[s]>_{\epsilon_i}$, this corresponds to the first part in the third term in the governning equation, and we omit the minus here.
            %> For this part, we note that the zero value of the interpolation function at the boundary point doesn't mean its derivative vanishes at the same point
            OP11(eidM, :) = OP11(eidM, :) - 0.5 * Js * mesh.InnerEdge.M * Dn1(eidM,:);
 
            %> $<[p_h],{\nabla_h s}>_{\epsilon_i}$, this corresponds to the second term in the governning equation, and we omit the minus here 
            OP11(:, eidM) = OP11(:, eidM) - 0.5 * ( Dn1(eidM,:) )' * ( Js * mesh.InnerEdge.M );
            
            EdgeMassMatrix(eidM, eidM) = Js * mesh.InnerEdge.M;
            QX11 = QX11 - ElementMassMatrix\( 0.5 * nx * EdgeMassMatrix );
            QY11 = QY11 - ElementMassMatrix\( 0.5 * ny * EdgeMassMatrix );
            
            % Actually Js of the adjacent cell should be used, here we
            % consider the mesh is conforming, so Js of the adjacent cell
            % is the same as the local value
            OP12 = zeros(Np);
            %> this part is used to indicate influence of local cell has on adjacent cell
            %> $(\tau^k [p_h],[s])_{\epsilon_i}$
            OP12(eidP, eidM) = - ( obj.Tau * ( Js * mesh.InnerEdge.M ) );
            %> $<[p_h],{\nabla_s s}>_{\epsilon_i}$ 
            OP12(:, eidM) = OP12(:, eidM) - 0.5*( (Dn2(eidP,:))'* ( Js * mesh.InnerEdge.M ));
            %> $<{\nabla_h p_h},[s]>_{\epsilon_i}$
            OP12(eidP, :) = OP12(eidP, :) + 0.5 * ( Js * mesh.InnerEdge.M ) * Dn1(eidM, :);
            
            EdgeMassMatrix = zeros(Np);
            AdjacentElementMassMatrix = diag(mesh.J(:,adjacentEle))*mesh.cell.M;
            EdgeMassMatrix(eidP, eidM) =  Js * mesh.InnerEdge.M;
            QX12 = - AdjacentElementMassMatrix\(0.5 * nx .* EdgeMassMatrix );
            QY12 = - AdjacentElementMassMatrix\(0.5 * ny .* EdgeMassMatrix );
            SecondOrderTerm( adjacentrows(:), cols(:) ) = AdjacentElementMassMatrix\( - OP12 );
            qx( adjacentrows(:), cols(:) ) = QX12;
            qy( adjacentrows(:), cols(:) ) = QY12;
    end
end
SecondOrderTerm( localrows(:), cols(:) ) = ElementMassMatrix\( - OP11 );
qx( localrows(:), cols(:) ) = QX11;
qy( localrows(:), cols(:) ) = QY11;

end