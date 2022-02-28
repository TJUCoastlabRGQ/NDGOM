classdef NdgVertLimiter3d < NdgVertLimiter
    %NDGVERTLIMITER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        BoundaryEdge
        BottomBoundaryEdge
        SurfaceBoundaryEdge
        Nv2d
        Nvc2d
    end
    
    properties
        % Inverse Vandmonde matrix corresponding to the first order, is used to
        % calculate the mode coefficient of the corrected vertex value
        InvVandVert
        % Vandmonde matrix corresponding to the first order, is used to
        % calculate the intepolation value at the intepolation point
        VandInterp
    end
    
    methods
        
        function obj = NdgVertLimiter3d( mesh )
            obj = obj@NdgVertLimiter( mesh );
            warning('off');
            obj.BoundaryEdge = struct(mesh.BoundaryEdge);
            obj.BottomBoundaryEdge = struct(mesh.BottomBoundaryEdge);
            obj.SurfaceBoundaryEdge = struct(mesh.SurfaceBoundaryEdge);
            warning('on');
            obj.Nv2d = mesh.mesh2d.Nv;
            obj.Nvc2d = obj.assembleVertexCellConnect2d( mesh );
            obj.InvVandVert = zeros(mesh.cell.Nv);
            obj.VandInterp = zeros(mesh.cell.Np, mesh.cell.Nv);
            if mesh.cell.type == enumStdCell.PrismTri
                ColIndex = [1, 2, mesh.cell.N + 2, mesh.cell.Nz * mesh.cell.Nph + 1,...
                    mesh.cell.Nz * mesh.cell.Nph + 2, mesh.cell.Nz * mesh.cell.Nph + mesh.cell.N + 2];
                RowIndex = [1, mesh.cell.N + 1, (mesh.cell.N + 1)*(mesh.cell.N + 2)/2,...
                    mesh.cell.Nz * mesh.cell.Nph + 1, mesh.cell.Nz * mesh.cell.Nph + mesh.cell.N + 1, ...
                    mesh.cell.Nz * mesh.cell.Nph + (mesh.cell.N + 1)*(mesh.cell.N + 2)/2];
            elseif mesh.cell.type == enumStdCell.PrismQuad
                ColIndex = [1, 2, mesh.cell.N + 2, mesh.cell.N + 3,...
                    (mesh.cell.N + 1)^2+1, (mesh.cell.N + 1)^2+2,...
                    (mesh.cell.N + 1)^2+mesh.cell.N + 2, (mesh.cell.N + 1)^2+mesh.cell.N + 3];
                RowIndex = [1, mesh.cell.N + 1, (mesh.cell.N + 1)*mesh.cell.N + 1, (mesh.cell.N + 1)^2,...
                    mesh.cell.Nz * mesh.cell.Nph + 1, mesh.cell.Nz * mesh.cell.Nph + mesh.cell.N + 1,...
                    mesh.cell.Nz * mesh.cell.Nph + (mesh.cell.N + 1)*mesh.cell.N + 1,...
                    mesh.cell.Nz * mesh.cell.Nph + (mesh.cell.N + 1)^2];
                
            end
            obj.InvVandVert = inv(mesh.cell.V(RowIndex, ColIndex));
            obj.VandInterp = mesh.cell.V(:, ColIndex);
        end
        
        function fphys = matLimit( obj, fphys, fieldId )
            % calculate the cell averages
            Nmesh = obj.Nmesh;
            cvar = cell( Nmesh, 1 );
            for m = 1:Nmesh
                cvar{m} = obj.meshUnion(m).GetMeshAverageValue(...
                    fphys{m}(:,:,fieldId) );
                fphys{m}(:,:,fieldId) = mxVertLimit3d( fphys{m}(:,:,fieldId), cvar{m},...
                    obj.meshUnion(m).EToE, inv(obj.meshUnion(m).cell.V1d),...
                    obj.meshUnion(m).cell.Npz, obj.meshUnion(m).cell.Nph, ...
                    obj.meshUnion(m).cell.V1d);
            end
        end
        
        function fphys = matLimitNew( obj, physClass, fphys )
            
%             fphys{1}(:,:,physClass.varFieldIndex) = ...
%                 mxVertLimit3dNew(fphys{1}, physClass.varFieldIndex, obj.meshUnion.cell.Np, ...
%                 obj.meshUnion.K, obj.meshUnion.J, obj.meshUnion.cell.wq, obj.meshUnion.cell.Vq, physClass.fext3d{1},...
%                 physClass.gra, physClass.hcrit, int8(physClass.meshUnion.BoundaryEdge.ftype), ...
%                 obj.BoundaryEdge, ...
%                 obj.BottomBoundaryEdge, obj.SurfaceBoundaryEdge, obj.Nv, obj.Nvc, obj.VToK,...
%                 obj.meshUnion.mesh2d.cell.Fmask, obj.meshUnion.mesh2d.cell.Np,...
%                 obj.meshUnion.cell.Nz, obj.meshUnion.cell.N, obj.meshUnion.mesh2d.cell.Nv, obj.meshUnion.LAV,...
%                 obj.meshUnion.EToV, obj.meshUnion.Nz, obj.Nv2d, obj.Nvc2d, obj.InvVandVert, obj.VandInterp);
            fphys{1}(:,:,physClass.varFieldIndex) = ...
                mxVertLimit3dWithLargeStencil(fphys{1}, physClass.varFieldIndex, obj.meshUnion.cell.Np, ...
                obj.meshUnion.K, obj.meshUnion.J, obj.meshUnion.cell.wq, obj.meshUnion.cell.Vq, ...
                obj.Nv, obj.Nvc, obj.VToK,...
                obj.meshUnion.mesh2d.cell.Fmask, obj.meshUnion.mesh2d.cell.Np,...
                obj.meshUnion.cell.Nz, obj.meshUnion.cell.N, obj.meshUnion.mesh2d.cell.Nv, obj.meshUnion.LAV,...
                obj.meshUnion.EToV, obj.meshUnion.Nz, obj.Nv2d, obj.Nvc2d, obj.InvVandVert, obj.VandInterp);
        end
    end
    
    methods( Access = protected )
        function Nvc = assembleVertexCellConnect2d( obj, mesh )
            Nv = mesh.mesh2d.Nv; % total number of vertex
            Nvc = zeros(Nv, 1); % number of cells connecting to each vertex
            for m = 1:obj.Nmesh
                for k = 1:mesh.mesh2d(m).K
                    v = mesh.mesh2d(m).EToV(:, k); % get the vertex index
                    Nvc(v) = Nvc(v) + 1;
                end
            end
        end
    end
    
end

