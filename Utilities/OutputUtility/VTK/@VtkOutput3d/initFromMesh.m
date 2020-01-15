function initFromMesh( obj, mesh2d, mesh3d )
    %initFromMesh - Description
    %
    % Syntax: output = initFromMesh(obj)
    %
    % Long description
    
    initFromMesh@VtkOutput2d(obj,mesh2d);
    [ obj.Np3d, Ncell, Ncon, obj.SEToV3d, ctype ] = InitSubConnect( mesh3d );

    obj.CellVertList3d = zeros( obj.Np3d, Ncell * mesh3d.K );

    sk = 1;
    for k = 1 : mesh3d.K
        nodesk = mesh3d.cell.Np * ( k - 1 );
        for m = 1 : Ncell
            obj.CellVertList3d( :, sk ) = nodesk + obj.SEToV3d(:, m) - 1; % transform to 0-offset
            sk = sk + 1;
        end
    end

    obj.Npoint3d = mesh3d.K * mesh3d.cell.Np;
    obj.Points3d = [ mesh3d.x(:), mesh3d.y(:), mesh3d.z(:) ]';
    obj.Ncell3d = mesh3d.K * Ncell;
    obj.Ncon3d = obj.Ncell3d * Ncon;
    obj.ctype3d = repmat( int8(ctype), obj.Ncell3d, 1 );
    
%     if ~isdir(obj.casename)
%         mkdir(obj.casename);
%     end
end

function [ Np, Ncell, Ncon, EToV, ctype ] = InitSubConnect( mesh )
    if ( mesh.cell.type == enumStdCell.PrismTri )
        [ Ncell, EToV ] = InitPrismTriConnect2d( mesh.cell.N, mesh.cell.Nz );
        Np = 6;
        Ncon = Np + 1;
        ctype = enumVtkCell.VTK_WEDGE;
    elseif ( mesh.cell.type == enumStdCell.PrismQuad )
        [ Ncell, EToV ] = InitPrismQuadConnect2d( mesh.cell.N, mesh.cell.Nz );
        Np = 8;
        Ncon = Np + 1;
        ctype = enumVtkCell.VTK_HEXAHEDRON;
    end
end