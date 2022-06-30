function [ mesh ] = makeSMSFileUMeshUnion2d( N, filename )
fid1 = fopen(filename, 'r');
if( fid1 < 0 )
    msgID = [mfilename, ':inputFileNameError'];
    msgtext = ['The input file name: ', filename, ' is incorrect'];
    ME = MException(msgID, msgtext);
    throw(ME);
end
% READ FIRST LINE
fgetl(fid1);

% READ POINT
temp = fscanf(fid1,'%d',2); Nv = temp(2); Ne = temp(1);
data = fscanf(fid1,'%d %f %f %f\n',[4,Nv]);
vx = data(2, :)';
vy = data(3, :)';

% READ CELL DATA

Nedge = 0;
Ntri = 0; EToVTri = cell(1); EToRTri = cell(1);
Nquad = 0;EToVQuad = cell(1); EToRQuad = cell(1);
for i = 1:Ne
    temp = fgetl(fid1);
    data = str2num(temp);
    if data(2) == 3
        Ntri = Ntri + 1;
        EToVTri{1}( Ntri, : ) = [data(3) data(4) data(5)];
        EToRTri{1}(Ntri) = 1;
    elseif data(2) == 4
        Nquad = Nquad + 1;
        EToVQuad{1}( Nquad, : ) = [data(3) data(4) data(5) data(6)];
        EToRQuad{1}(Ntri) = 1;
    else
        msgID = [mfilename, ':CellTypeError'];
        msgtext = ['The cell type: ', data(2), ' is unknown'];
        ME = MException(msgID, msgtext);
        throw(ME);
    end
end

if Ntri > 0
    stdTri = StdTri(N);
    triMesh = NdgMesh2d(stdTri, Nv, vx, vy, Ntri, (EToVTri{1})', EToRTri{1});
end

if Nquad > 0
    stdQuad = StdQuad(N);
    quadMesh = NdgMesh2d(stdQuad, Nv, vx, vy, Nquad, EToVQuad{1}, EToRQuad{1});
end

BCToV = cell(1); IBCToV = 1;
str = fgetl(fid1);
k = strfind(str,'Number of open boundaries');
if k
    tempstr=isstrprop(str,'digit');
    tempstr=str(tempstr);
    Nopen = str2num(tempstr);
    fgetl(fid1); % The rest part of the line include TNOP
    for boundary = 1 : Nopen
        TempNOP =  fscanf(fid1,'%d',1);
        fgetl(fid1); % The rest part of the line include NOP
        IOOBN = fscanf(fid1,'%d',TempNOP); % Index of open boundary nodes
        for i = 1:numel(IOOBN) - 1
%             BCToV{1}( IBCToV, : ) = [ IOOBN(i), IOOBN(i+1), 6]; % 6 stands for clamped depth condition
            BCToV{1}( :, IBCToV ) = [ IOOBN(i), IOOBN(i+1), 6];
            IBCToV = IBCToV + 1;
        end
    end
end

% fgetl(fid1); % Maybe a blankspace
str = fgetl(fid1);
k = strfind(str,'Number of land boundaries');
if k
    tempstr=isstrprop(str,'digit');
    tempstr=str(tempstr);
    NOW = str2num(tempstr); % Number of slip wall
    fgetl(fid1); % The rest part of the line include TNOP
    for boundary = 1 : NOW
        TempNOP =  fscanf(fid1,'%d',1);
        fgetl(fid1); % The rest part of the line include NOP
        IOSWN = fscanf(fid1,'%d',TempNOP); % Index of slip wall boundary nodes
        for i = 1:numel(IOSWN) - 1
%             BCToV{1}( IBCToV, : ) = [ IOSWN(i), IOSWN(i+1), 2]; % 6 stands for slip condition
            BCToV{1}( :, IBCToV ) = [ IOSWN(i), IOSWN(i+1), 2]; % 6 stands for slip condition
            IBCToV = IBCToV + 1;
        end
    end
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
    mesh.BoundaryEdge = NdgHaloEdge2d( mesh, mesh.ind, BCToV{1} );
elseif (Nquad > 0) && (Ntri == 0)
    mesh = quadMesh;
    mesh.ConnectMeshUnion( 1, mesh);
    mesh.InnerEdge = NdgInnerEdge2d( mesh, 1 );
    mesh.BoundaryEdge = NdgHaloEdge2d( mesh, mesh.ind, BCToV );
end

end

