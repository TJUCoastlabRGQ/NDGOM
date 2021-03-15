function [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d( mesh2d, mesh3d , str)
%%periodic boundary condition for this case
if (mesh3d.cell.type == enumStdCell.PrismQuad )
    if strcmp(str,'West-East')
        LOrder2d = find(all(mesh2d.BoundaryEdge.xb == min(min(mesh2d.BoundaryEdge.xb))));
        LOrder3d = find(all(mesh3d.BoundaryEdge.xb == min(min(mesh3d.BoundaryEdge.xb))));
        ROrder2d = find(all(mesh2d.BoundaryEdge.xb == max(max(mesh2d.BoundaryEdge.xb))));
        ROrder3d = find(all(mesh3d.BoundaryEdge.xb == max(max(mesh3d.BoundaryEdge.xb))));
    elseif strcmp(str,'South-North')%Rotate clockwise
        LOrder2d = find(all(mesh2d.BoundaryEdge.yb == min(min(mesh2d.BoundaryEdge.yb))));
        LOrder3d = find(all(mesh3d.BoundaryEdge.yb == min(min(mesh3d.BoundaryEdge.yb))));
        ROrder2d = find(all(mesh2d.BoundaryEdge.yb == max(max(mesh2d.BoundaryEdge.yb))));
        ROrder3d = find(all(mesh3d.BoundaryEdge.yb == max(max(mesh3d.BoundaryEdge.yb))));
    elseif strcmp(str,'Surface-Bottom')
        %This is to be finished later
    end
    %Ne part, this is the same for the west-east and the south-north boundary
    
    mesh2d.BoundaryEdge.Ne = mesh2d.BoundaryEdge.Ne - numel(LOrder2d) - numel(ROrder2d);
    mesh3d.BoundaryEdge.Ne = mesh3d.BoundaryEdge.Ne - numel(LOrder3d) - numel(ROrder3d);
    mesh2d.InnerEdge.Ne = mesh2d.InnerEdge.Ne + numel(LOrder2d);
    mesh3d.InnerEdge.Ne = mesh3d.InnerEdge.Ne + numel(LOrder3d);
    
    %ftype part, this is the same for the west-east and the south-north boundary
    mesh2d.BoundaryEdge.ftype([LOrder2d,ROrder2d])=[];
    mesh3d.BoundaryEdge.ftype([LOrder3d,ROrder3d])=[];
    
    %EToE, FToE and FToF part, this part is not the same
    InnerEdgeFToE2d = zeros(2, mesh2d.InnerEdge.Ne);
    InnerEdgeFToE3d = zeros(2, mesh3d.InnerEdge.Ne);
    InnerEdgeFToF2d = zeros(2, mesh2d.InnerEdge.Ne);
    InnerEdgeFToF3d = zeros(2, mesh3d.InnerEdge.Ne);
    
    if strcmp(str,'West-East')
        LYCor = sort(mesh2d.BoundaryEdge.yb(:,LOrder2d));
        RYCor = sort(mesh2d.BoundaryEdge.yb(:,ROrder2d));
    elseif strcmp(str,'South-North')
        LYCor = sort(mesh2d.BoundaryEdge.xb(:,LOrder2d));
        RYCor = sort(mesh2d.BoundaryEdge.xb(:,ROrder2d));
    elseif strcmp(str,'Surface-Bottom')
        %This is to be finished later
    end
    for i = 1:size(LYCor,2)
        order = find(all(bsxfun(@eq,RYCor,LYCor(:,i)),1));
        ele1 = mesh2d.BoundaryEdge.FToE(1,LOrder2d(i));
        face1 = mesh2d.BoundaryEdge.FToF(1,LOrder2d(i));
        ele2 = mesh2d.BoundaryEdge.FToE(1,ROrder2d(order));
        face2 = mesh2d.BoundaryEdge.FToF(1,ROrder2d(order));
        mesh2d.EToE((ele2-1)*mesh2d.cell.Nface+face2) = ele1;
        mesh2d.EToE((ele1-1)*mesh2d.cell.Nface+face1) = ele2;
        InnerEdgeFToE2d(1,i) = ele2;
        InnerEdgeFToF2d(1,i) = face2;
        InnerEdgeFToE2d(2,i) = ele1;
        InnerEdgeFToF2d(2,i) = face1;
        for j = 1:mesh3d.Nz
            InnerEdgeFToE3d(1,(i-1)*mesh3d.Nz + j) = (ele2-1)*mesh3d.Nz + j;
            InnerEdgeFToE3d(2,(i-1)*mesh3d.Nz + j) = (ele1-1)*mesh3d.Nz + j;
            InnerEdgeFToF3d(1,(i-1)*mesh3d.Nz + j) = face2;
            InnerEdgeFToF3d(2,(i-1)*mesh3d.Nz + j) = face1;
            mesh3d.EToE(((ele2-1)*mesh3d.Nz + j-1)*mesh3d.cell.Nface+face2) = (ele1-1)*mesh3d.Nz + j;
            mesh3d.EToE(((ele1-1)*mesh3d.Nz + j-1)*mesh3d.cell.Nface+face1) = (ele2-1)*mesh3d.Nz + j;
        end
    end
    InnerEdgeFToE3d(:,numel(LOrder3d)+1:end) = mesh3d.InnerEdge.FToE;
    InnerEdgeFToE2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.FToE;
    InnerEdgeFToF3d(:,numel(LOrder3d)+1:end) = mesh3d.InnerEdge.FToF;
    InnerEdgeFToF2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.FToF;
    mesh3d.InnerEdge.FToE = InnerEdgeFToE3d;
    mesh2d.InnerEdge.FToE = InnerEdgeFToE2d;
    mesh3d.InnerEdge.FToF = InnerEdgeFToF3d;
    mesh2d.InnerEdge.FToF = InnerEdgeFToF2d;
    mesh2d.BoundaryEdge.FToE(:,[LOrder2d, ROrder2d]) = [];
    mesh3d.BoundaryEdge.FToE(:,[LOrder3d, ROrder3d]) = [];
    mesh2d.BoundaryEdge.FToF(:,[LOrder2d, ROrder2d]) = [];
    mesh3d.BoundaryEdge.FToF(:,[LOrder3d, ROrder3d]) = [];
    %nx, ny and nz part
    InnerEdgenx3d = zeros(size(mesh3d.InnerEdge.nx,1),mesh3d.InnerEdge.Ne);
    InnerEdgeny3d = zeros(size(mesh3d.InnerEdge.ny,1),mesh3d.InnerEdge.Ne);
    InnerEdgenz3d = zeros(size(mesh3d.InnerEdge.nz,1),mesh3d.InnerEdge.Ne);
    InnerEdgenx3d(:,1:numel(ROrder3d)) = mesh3d.BoundaryEdge.nx(:,ROrder3d);
    InnerEdgenx3d(:,numel(ROrder3d)+1:end) = mesh3d.InnerEdge.nx;
    InnerEdgeny3d(:,1:numel(ROrder3d)) = mesh3d.BoundaryEdge.ny(:,ROrder3d);
    InnerEdgeny3d(:,numel(ROrder3d)+1:end) = mesh3d.InnerEdge.ny;
    InnerEdgenz3d(:,1:numel(ROrder3d)) = mesh3d.BoundaryEdge.nz(:,ROrder3d);
    InnerEdgenz3d(:,numel(ROrder3d)+1:end) = mesh3d.InnerEdge.nz;
   
    InnerEdgenx2d = zeros(size(mesh2d.InnerEdge.nx,1),mesh2d.InnerEdge.Ne);
    InnerEdgeny2d = zeros(size(mesh2d.InnerEdge.ny,1),mesh2d.InnerEdge.Ne);
    InnerEdgenx2d(:,1:numel(ROrder2d)) = mesh2d.BoundaryEdge.nx(:,ROrder2d);
    InnerEdgenx2d(:,numel(ROrder2d)+1:end) = mesh2d.InnerEdge.nx;
    InnerEdgeny2d(:,1:numel(ROrder2d)) = mesh2d.BoundaryEdge.ny(:,ROrder2d);
    InnerEdgeny2d(:,numel(ROrder2d)+1:end) = mesh2d.InnerEdge.ny;
    
    mesh2d.InnerEdge.nx = InnerEdgenx2d; mesh2d.InnerEdge.ny = InnerEdgeny2d;
    mesh3d.InnerEdge.nx = InnerEdgenx3d; mesh3d.InnerEdge.ny = InnerEdgeny3d;
    mesh3d.InnerEdge.nz = InnerEdgenz3d;
    mesh2d.BoundaryEdge.nx(:,[LOrder2d, ROrder2d]) = [];
    mesh2d.BoundaryEdge.ny(:,[LOrder2d, ROrder2d]) = [];
    mesh3d.BoundaryEdge.nx(:,[LOrder3d, ROrder3d]) = [];
    mesh3d.BoundaryEdge.ny(:,[LOrder3d, ROrder3d]) = [];
    mesh3d.BoundaryEdge.nz(:,[LOrder3d, ROrder3d]) = [];
    %xb, yb and zb part
    mesh2d.BoundaryEdge.xb(:,[LOrder2d, ROrder2d]) = [];
    mesh2d.BoundaryEdge.yb(:,[LOrder2d, ROrder2d]) = [];
    mesh3d.BoundaryEdge.xb(:,[LOrder3d, ROrder3d]) = [];
    mesh3d.BoundaryEdge.yb(:,[LOrder3d, ROrder3d]) = [];
    mesh3d.BoundaryEdge.zb(:,[LOrder3d, ROrder3d]) = [];
    %LAV part(for 2d mesh only)
    InnerEdgeLAV2d = zeros(size(mesh2d.InnerEdge.LAV,1),mesh2d.InnerEdge.Ne);
    InnerEdgeLAV2d(:,1:numel(ROrder2d)) = mesh2d.BoundaryEdge.LAV(:,ROrder2d);
    InnerEdgeLAV2d(:,numel(ROrder2d)+1:end) = mesh2d.InnerEdge.LAV;
    mesh2d.InnerEdge.LAV = InnerEdgeLAV2d;
    mesh2d.BoundaryEdge.LAV(:,[LOrder2d, ROrder2d]) = [];
    
    InnerEdgeLAV3d = zeros(size(mesh3d.InnerEdge.LAV,1),mesh3d.InnerEdge.Ne);
    InnerEdgeLAV3d(:,1:numel(ROrder3d)) = mesh3d.BoundaryEdge.LAV(:,ROrder3d);
    InnerEdgeLAV3d(:,numel(ROrder3d)+1:end) = mesh3d.InnerEdge.LAV;    
    mesh3d.InnerEdge.LAV = InnerEdgeLAV3d;
    mesh3d.BoundaryEdge.LAV(:,[LOrder3d, ROrder3d]) = [];
    
    %FToN1 and FToN2 part, this is the same for the west-east and the south-north boundary
    InnerEdgeFToN13d = zeros(size(mesh3d.InnerEdge.FToN1,1),mesh3d.InnerEdge.Ne);
    InnerEdgeFToN23d = zeros(size(mesh3d.InnerEdge.FToN2,1),mesh3d.InnerEdge.Ne);
    InnerEdgeFToN12d = zeros(size(mesh2d.InnerEdge.FToN1,1),mesh2d.InnerEdge.Ne);
    InnerEdgeFToN22d = zeros(size(mesh2d.InnerEdge.FToN2,1),mesh2d.InnerEdge.Ne);
    InnerEdgeFToN12d(:,1:numel(ROrder2d)) = mesh2d.BoundaryEdge.FToN1(:,ROrder2d);
    InnerEdgeFToN22d(:,1:numel(ROrder2d)) = flip(mesh2d.BoundaryEdge.FToN2(:,LOrder2d));
    InnerEdgeFToN13d(:,1:numel(ROrder3d)) = mesh3d.BoundaryEdge.FToN1(:,ROrder3d);
    InnerEdgeFToN23d(:,1:numel(ROrder3d)) = mesh3d.BoundaryEdge.FToN2(:,LOrder3d);
    InnerEdgeFToN12d(:,numel(ROrder2d)+1:end) = mesh2d.InnerEdge.FToN1;
    InnerEdgeFToN22d(:,numel(ROrder2d)+1:end) = mesh2d.InnerEdge.FToN2;
    InnerEdgeFToN13d(:,numel(ROrder3d)+1:end) = mesh3d.InnerEdge.FToN1;
    InnerEdgeFToN23d(:,numel(ROrder3d)+1:end) = mesh3d.InnerEdge.FToN2;
    mesh2d.InnerEdge.FToN1 = InnerEdgeFToN12d;
    mesh2d.InnerEdge.FToN2 = InnerEdgeFToN22d;
    mesh3d.InnerEdge.FToN1 = InnerEdgeFToN13d;
    mesh3d.InnerEdge.FToN2 = InnerEdgeFToN23d;
    mesh2d.BoundaryEdge.FToN1(:,[LOrder2d, ROrder2d]) = [];
    mesh2d.BoundaryEdge.FToN2(:,[LOrder2d, ROrder2d]) = [];
    mesh3d.BoundaryEdge.FToN1(:,[LOrder3d, ROrder3d]) = [];
    mesh3d.BoundaryEdge.FToN2(:,[LOrder3d, ROrder3d]) = [];
    %GFToN1 and GFToN2 part
    mesh3d.InnerEdge.GFToN1 = bsxfun(@plus, ( mesh3d.InnerEdge.FToE(1,:) - 1 ) * mesh3d.cell.Np , mesh3d.InnerEdge.FToN1);
    mesh3d.InnerEdge.GFToN2 = bsxfun(@plus, ( mesh3d.InnerEdge.FToE(2,:) - 1 ) * mesh3d.cell.Np , mesh3d.InnerEdge.FToN2);
    mesh3d.BoundaryEdge.GFToN1 = bsxfun(@plus, ( mesh3d.BoundaryEdge.FToE(1,:) - 1 ) * mesh3d.cell.Np , mesh3d.BoundaryEdge.FToN1);
    mesh3d.BoundaryEdge.GFToN2 = bsxfun(@plus, ( mesh3d.BoundaryEdge.FToE(2,:) - 1 ) * mesh3d.cell.Np , mesh3d.BoundaryEdge.FToN2);
    mesh3d.InnerEdge.Jz = mesh3d.Jz(mesh3d.InnerEdge.GFToN1);
    mesh3d.BoundaryEdge.Jz = mesh3d.Jz(mesh3d.BoundaryEdge.GFToN1);
    %Js part
    InnerEdgeJs3d = zeros(size(mesh3d.InnerEdge.Js,1),mesh3d.InnerEdge.Ne);
    InnerEdgeJs2d = zeros(size(mesh2d.InnerEdge.Js,1),mesh2d.InnerEdge.Ne);
    InnerEdgeJs2d(:,1:numel(ROrder2d)) = mesh2d.BoundaryEdge.Js(:,ROrder2d);
    InnerEdgeJs3d(:,1:numel(ROrder3d)) = mesh3d.BoundaryEdge.Js(:,ROrder3d);
    InnerEdgeJs2d(:,numel(ROrder2d)+1:end) = mesh2d.InnerEdge.Js;
    InnerEdgeJs3d(:,numel(ROrder3d)+1:end) = mesh3d.InnerEdge.Js;
    mesh2d.InnerEdge.Js = InnerEdgeJs2d;
    mesh3d.InnerEdge.Js = InnerEdgeJs3d;
    mesh2d.BoundaryEdge.Js(:,[LOrder2d, ROrder2d]) = [];
    mesh3d.BoundaryEdge.Js(:,[LOrder3d, ROrder3d]) = [];
else
    msgID = 'ImposePeriodicBoundaryCondition:CellTypeError';
    msgtext = 'Periodic boundary condition is only applicable for quad and prismquad.';
    ME = MException(msgID, msgtext);
    throw(ME);    
end
end