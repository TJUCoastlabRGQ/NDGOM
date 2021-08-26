function [ mesh ] = makeGmshFileUMeshUnion2d( N, Elementfilename, Boundaryfilename )
%> This file is updated on 2021/07/09 by Guoquan Ran to accomadate the storement manner
%> of date by Gmesh since the boundary condition cannot be included in one single file, 
%> i.e. "Elementfilename". The mesh is created with GMSH-4.8.4, while the data is stored
%> in version 2(File -> Export -> Version 2 ASCLL). Difference between the "Elementfilename" 
%> and the "Boundaryfilename" is that the choice, Save all elements, is included for the
%> former file but not included for the latter file.
%> NOTE: The format of the boundary edge date is as follows:
%> elm-number elm-type number-of-tags physical-tag element-tag node-number-list
fid1 = fopen(Elementfilename, 'r');
fid2 = fopen(Boundaryfilename, 'r');
if( fid1 < 0 || fid2 < 0 )
    msgID = [mfilename, ':inputFileNameError'];
    msgtext = ['The input file name: ', Elementfilename, 'or :', Boundaryfilename, ' is incorrect'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
% Jump the mesh format part for both file
for i=1:4
   fgetl(fid1);
   fgetl(fid2);
end

Nv = fscanf(fid1,'%d\n',1);
BNv = fscanf(fid2,'%d\n',1);
%Jump the node part for boundary file
for i=1:BNv
    fgetl(fid2);
end
data = fscanf(fid1,'%d %f %f %f\n',[4,Nv]);
vx = data(2, :)';
vy = data(3, :)';
%Jump the tag part and enter the element part
for i=1:2
   fgetl(fid1);
   fgetl(fid2);
end
Ne = fscanf(fid1,'%d\n',1);
EleNedge = 0;
Nedge = fscanf(fid2,'%d\n',1);
Ntri = 0; 
Nquad = 0;
Npoint = 0;
for i = 1:Ne
    temp = fgetl(fid1);
    data = str2num(temp);
    switch data(2)
        case 1
            EleNedge = EleNedge + 1;
        case 2
            Ntri = Ntri + 1;
        case 3
            Nquad = Nquad + 1;
        case 15
            Npoint = Npoint + 1;
    end
end

fseek(fid1,0,-1);
for i=1:(Nv+8+Npoint + EleNedge)
    fgetl(fid1);
end

data = fscanf(fid2,'%d %d %d %d %d %d %d\n',[7,Nedge]);
BCToV = [data(6, :); data(7, :); data(4, :)];

if Ntri > 0
    data = fscanf(fid1,'%d %d %d %d %d %d %d %d \n',[8, Ntri]);
    EToVTri = [data(6, :); data(7, :); data(8, :)];
    EToRTri = data(4, :);
    stdTri = StdTri(N);
    triMesh = NdgMesh2d(stdTri, Nv, vx, vy, Ntri, EToVTri, EToRTri);
end

if Nquad > 0
    data = fscanf(fid1,'%d %d %d %d %d %d %d %d %d\n',[9, Nquad]);
    EToVQuad = [data(6, :); data(7, :); data(8, :); data(9, :)];
    EToRQuad = data(4, :);
    stdQuad = StdQuad(N);
    quadMesh = NdgMesh2d(stdQuad, Nv, vx, vy, Nquad, EToVQuad, EToRQuad);
end

if (Ntri > 0) && (Nquad > 0)
    mesh = [triMesh, quadMesh];
    mesh(1).ConnectMeshUnion( 1, mesh);
    mesh(1).InnerEdge = NdgInnerEdge2d( mesh, mesh(1).ind );
    mesh(1).BoundaryEdge = NdgHaloEdge2d( mesh, mesh(1).ind, BCToV );
    
    mesh(2).ConnectMeshUnion( 2, mesh);
    mesh(2).InnerEdge = NdgInnerEdge2d( mesh, mesh(2).ind );
    mesh(2).BoundaryEdge = NdgHaloEdge2d( mesh, mesh(2).ind, BCToV );
elseif (Ntri > 0) && (Nquad == 0)
    mesh = triMesh;
    mesh.ConnectMeshUnion( 1, mesh);
    mesh.InnerEdge = NdgInnerEdge2d( mesh, 1 );
    mesh.BoundaryEdge = NdgHaloEdge2d( mesh, mesh.ind, BCToV );
elseif (Nquad > 0) && (Ntri == 0)
    mesh = quadMesh;
    mesh.ConnectMeshUnion( 1, mesh);
    mesh.InnerEdge = NdgInnerEdge2d( mesh, 1 );
    mesh.BoundaryEdge = NdgHaloEdge2d( mesh, mesh.ind, BCToV );
end

fclose(fid1);
fclose(fid2);
end

