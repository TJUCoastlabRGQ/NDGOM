function [ mesh2d ] = ImposePeriodicBoundaryCondition2d( mesh2d, str)
%we note that FTOF property is not corrected here, and this may be needed
%for non-hydrostatic nodel
%%periodic boundary condition for this case
if strcmp(str,'West-East')
    LOrder2d = find(all(mesh2d.BoundaryEdge.xb == min(min(mesh2d.BoundaryEdge.xb))));
    ROrder2d = find(all(mesh2d.BoundaryEdge.xb == max(max(mesh2d.BoundaryEdge.xb))));
elseif strcmp(str,'South-North')
    LOrder2d = find(all(mesh2d.BoundaryEdge.yb == min(min(mesh2d.BoundaryEdge.yb))));
    ROrder2d = find(all(mesh2d.BoundaryEdge.yb == max(max(mesh2d.BoundaryEdge.yb))));
end
%Ne part, this is the same for the west-east and the south-north boundary

mesh2d.BoundaryEdge.Ne = mesh2d.BoundaryEdge.Ne - numel(LOrder2d) - numel(ROrder2d);
mesh2d.InnerEdge.Ne = mesh2d.InnerEdge.Ne + numel(LOrder2d);

%ftype part, this is the same for the west-east and the south-north boundary
mesh2d.BoundaryEdge.ftype([LOrder2d,ROrder2d])=[];

%FToE part and FToF part, this part is not the same
InnerEdgeFToE2d = zeros(2, mesh2d.InnerEdge.Ne);
InnerEdgeFToF2d = zeros(2, mesh2d.InnerEdge.Ne);

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
    InnerEdgeFToE2d(1,i) = ele1;
    InnerEdgeFToF2d(1,i) = face1;
    InnerEdgeFToE2d(2,i) = ele2;
    InnerEdgeFToF2d(2,i) = face2;
end
InnerEdgeFToE2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.FToE;
InnerEdgeFToF2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.FToF;
mesh2d.InnerEdge.FToE = InnerEdgeFToE2d;
mesh2d.InnerEdge.FToF = InnerEdgeFToF2d;
mesh2d.BoundaryEdge.FToE(:,[LOrder2d, ROrder2d]) = [];
mesh2d.BoundaryEdge.FToF(:,[LOrder2d, ROrder2d]) = [];
%nx, ny part

InnerEdgenx2d = zeros(size(mesh2d.InnerEdge.nx,1),mesh2d.InnerEdge.Ne);
InnerEdgeny2d = zeros(size(mesh2d.InnerEdge.ny,1),mesh2d.InnerEdge.Ne);
InnerEdgenx2d(:,1:numel(LOrder2d)) = mesh2d.BoundaryEdge.nx(:,LOrder2d);
InnerEdgenx2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.nx;
InnerEdgeny2d(:,1:numel(LOrder2d)) = mesh2d.BoundaryEdge.ny(:,LOrder2d);
InnerEdgeny2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.ny;

mesh2d.InnerEdge.nx = InnerEdgenx2d; mesh2d.InnerEdge.ny = InnerEdgeny2d;
mesh2d.BoundaryEdge.nx(:,[LOrder2d, ROrder2d]) = [];
mesh2d.BoundaryEdge.ny(:,[LOrder2d, ROrder2d]) = [];

%xb, yb and zb part
mesh2d.BoundaryEdge.xb(:,[LOrder2d, ROrder2d]) = [];
mesh2d.BoundaryEdge.yb(:,[LOrder2d, ROrder2d]) = [];

%LAV part(for 2d mesh only)
InnerEdgeLAV2d = zeros(size(mesh2d.InnerEdge.LAV,1),mesh2d.InnerEdge.Ne);
InnerEdgeLAV2d(:,1:numel(LOrder2d)) = mesh2d.BoundaryEdge.LAV(:,LOrder2d);
InnerEdgeLAV2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.LAV;
mesh2d.InnerEdge.LAV = InnerEdgeLAV2d;
mesh2d.BoundaryEdge.LAV(:,[LOrder2d, ROrder2d]) = [];


%FToN1 and FToN2 part, this is the same for the west-east and the south-north boundary
InnerEdgeFToN12d = zeros(size(mesh2d.InnerEdge.FToN1,1),mesh2d.InnerEdge.Ne);
InnerEdgeFToN22d = zeros(size(mesh2d.InnerEdge.FToN2,1),mesh2d.InnerEdge.Ne);
InnerEdgeFToN12d(:,1:numel(LOrder2d)) = mesh2d.BoundaryEdge.FToN1(:,LOrder2d);
InnerEdgeFToN22d(:,1:numel(ROrder2d)) = flip(mesh2d.BoundaryEdge.FToN2(:,ROrder2d));
InnerEdgeFToN12d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.FToN1;
InnerEdgeFToN22d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.FToN2;
mesh2d.InnerEdge.FToN1 = InnerEdgeFToN12d;
mesh2d.InnerEdge.FToN2 = InnerEdgeFToN22d;
mesh2d.BoundaryEdge.FToN1(:,[LOrder2d, ROrder2d]) = [];
mesh2d.BoundaryEdge.FToN2(:,[LOrder2d, ROrder2d]) = [];

%Js part
InnerEdgeJs2d = zeros(size(mesh2d.InnerEdge.Js,1),mesh2d.InnerEdge.Ne);
InnerEdgeJs2d(:,1:numel(LOrder2d)) = mesh2d.BoundaryEdge.Js(:,LOrder2d);
InnerEdgeJs2d(:,numel(LOrder2d)+1:end) = mesh2d.InnerEdge.Js;
mesh2d.InnerEdge.Js = InnerEdgeJs2d;
mesh2d.BoundaryEdge.Js(:,[LOrder2d, ROrder2d]) = [];
end